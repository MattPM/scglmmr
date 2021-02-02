#' GroupContrastGLMMtidy - fit a random intercept model of the difference in fold change between groups and fold change across both groups and baseline controlling for any number of covariates.
#'
#' @param tidydf tidy data must contain column 'subjectid' and any names covariates (specified with )
#' @param response_variable the column containing the variables to be modeled in tidy (long) format -- e.g. 'gene' with the gene's values in the column in response_value
#' @param response_value the value for the measurement response_variable -- should probably be in log space.
#' @param fixed_effects the fixed effects to specify in this case as covariates e.g. c('gender', 'ethnicity'), these are automatically combined into formula with your_fixed_effects + 0 + group_id + (1|subjectid)
#' @param lmer_formula lmer model formula, does NOT need to be specified. If NULL and if fixed_effects are NULL, automatically creates a donor random intercept model of 0 + group_id + (1|subjectid)
#' @param plotdatqc if TRUE , saves a qq plot of the variable, random effects, variance covariance matrix, and em means vs data means boxplot fo the directory specified with figpath
#' @param figpath path to save figures e.g.  file.path("mypath/")
#' @return
#' @import ggplot2
#' @import lme4
#' @importFrom stats as.formula
#' @importFrom rlang sym
#' @importFrom emmeans emmeans
#' @importFrom lattice qqmath levelplot dotplot
#' @importFrom gridExtra grid.arrange
#' @importFrom viridis plasma
#' @importFrom ggpubr ggpaired rremove
#' @importFrom egg ggarrange
#' @importFrom broom tidy
#' @importFrom dplyr select mutate bind_rows if_else filter
#' @importFrom tidyr spread everything
#' @export
#'
#' @examples
#' \donttest{
#'
#'  # subjectid and group_id need to be specified, other variables can be added as covariates by specifying a vector of column names in tidydf to `fixed_effects`
#'  tidydf = mydata %>%
#'    dplyr::rename('subjectid'= subject) %>%
#'    mutate(group_id = factor(group_id, levels = c(0_0, 0_1, 1_0, 1_1))) %>%
#'    mutate(log10protvalue = log10(expression + 1)) %>%
#'
#'
#'  # tidydf should have a column with the response variable (i.e. variable on the formua eft hand side) in long form, i.e. a variable called "genes" or "proteins" or "cytokines" etc.
#'  # the values for the response variable are specified by the `response_value` for example, the log transfomred gene expression
#'
#'   mmres = GroupContrastGLMMtidy(tidydf = tidydf, response_variable = 'protein', response_value = 'log10protvalue',
#'    fixed_effects = c('gender', 'ethnicity'), lmer_formula = NULL, subject_colum = 'subjectid', plotdatqc  = TRUE, figpath = '~' )
#'   }
#'
GroupContrastGLMMtidy = function(tidydf, response_variable, response_value,
                                 fixed_effects = NULL, lmer_formula = NULL, plotdatqc = TRUE, figpath){

  # specify custom contrasts difference in treatment between groups, treatment effect across groups, per-treatment difference between groups.
  c00 = c(1,0,0,0) ; c01 = c(0,1,0,0) ; c10 = c(0,0,1,0) ; c11 = c(0,0,0,1)
  contrast_list = list(
    "time1vs0_group2vs1" = (c11 - c10) - (c01 - c00),
    "time1vs0" = (c11 + c01) / 2 - (c10 + c00) / 2,
    "time0_group2vs1" = c10 - c00
    )

  # specify model formula
  if(!is.null(lmer_formula)){
    f1 = lmer_formula
    } else if(!is.null(fixed_effects)) {
      fixef = paste0(fixed_effects, collapse = " + ")
      f1 = paste(response_value, "~ 0 +", fixef, '+ group_id + (1|subjectid)')
    } else {
      f1 = paste(response_value, "~ 0 +", 'group_id + (1|subjectid)')
    }
  f_store = f1
  f1 = stats::as.formula(f1)

  # checks on minimal required metadata and formatting of group levels
  print('model specified:'); print(f1)
  stopifnot(is.factor(metadata$group_id))
  cat("check levels, testing: \n",
      "1: difference in treatmant effect between groups: ",
      levels(metadata$group_id)[3:4], " fold change  vs: ", levels(metadata$group_id)[1:2], "fold change  \n",
      "2: pre treatment baseline difference between groups: ", levels(metadata$group_id)[c(1,3)], " fold change \n",
      "3: treatment effect across both groups combined \n",
      "a prori contrasts across `group_id` variable: "
  )
  print(contrast_list)


  # init store
  res_list = list()
  varsall = unique(tidydf[[response_variable]])
  # fit model for each protein
  for (i in 1:length(varsall)) {

    # fit model separately for each modeled  response variable
    gvar = rlang::sym('group_id')
    m1 = lme4::lmer(formula = f1,  data = tidydf[tidydf[[response_variable]] == varsall[i], ], REML = TRUE)
    emm1 = tryCatch(emmeans::emmeans(m1, specs = ~ {{ group_id }}), error = function(e)  return(NA))

    # error checking
    if( suppressWarnings(is.na(m1))) {
      print(paste0(" could not fit  model for ", varsall[i]))
      res_list[[i]] = NA
    }
    if(is.na(emm1)){
      print(paste0(varsall[i],  " check cov matrix could not calculate marginal means  "))
      res_list[[i]] = NA
    } else{

      # format model results
      marginal_mean_df = broom::tidy(emm1)[ ,1:2] %>% tidyr::spread('group_id','estimate')
      names(marginal_mean_df) = paste0(names(marginal_mean_df), "_marginal_mean")

      contrast_fit = emmeans::contrast(emm1, method = contrast_list)
      contrastdf = broom::tidy(contrast_fit)
      contrast_names = names(contrast_list)
      contrast_colnames = sapply(contrast_names, function(x) paste0(colnames(contrastdf), x)) %>%
        as.data.frame() %>%
        unlist(use.names = F) %>%
        as.character()
      # single line for feature i ; to be rbound
      contrastdf = cbind(contrastdf[1, ], contrastdf[2, ], contrastdf[3, ])
      names(contrastdf) = contrast_colnames

      # flag singular fits
      singular_fit = m1@optinfo$conv$lme4$messages
      contrastdf = contrastdf %>% dplyr::mutate(singular_fit = dplyr::if_else(is.null(singular_fit), true = 0, false = 1))

      # merge models
      res_list[[i]] =
        cbind(contrastdf, marginal_mean_df) %>%
        dplyr::mutate(variable = varsall[i]) %>%
        dplyr::mutate(formula = f_store) %>%
        dplyr::select(variable, tidyr::everything())

      # plot model quality control and data vs model view save to figpath
      if(isTRUE(plotdatqc)){
        # plot model quality control
        pdf(paste0(figpath,"qc ", varsall[i],".pdf"), width = 9, height = 9)
        p1 = lattice::qqmath(m1)
        p4 = lattice::levelplot(as.matrix(lme4::vcov.merMod(m1)),
                                main = "variance covariance matrix",
                                scales=list(y=list(cex=.5),x=list(rot=90, cex=.5)),
                                col.regions = viridis::plasma(100) )
        p2 = lattice::dotplot(lme4::ranef(m1))
        p3 = lattice::qqmath(lme4::ranef(m1))
        gridExtra::grid.arrange(p1,p2$subjectid,p3$subjectid,p4)
        dev.off()

        # plot model marginal means and boxplot of data means
        p1 =
          plot(emm1) +
          coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), axis.title.x = element_blank()) +
          xlab("model estimated marginal means") +
          scale_x_continuous(position = "left") +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

        p2 =
          ggpubr::ggpaired(tidydf[tidydf[[response_variable]] == varsall[i], ],x = 'group_id', y = response_value, fill = 'group_id') +
          scale_fill_manual(values = c("dodgerblue", "midnightblue", "red", "firebrick")) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90)) +
          ggtitle(varsall[i]) +
          ylab("modeled value") +
          xlab("group_id") +
          ggpubr::rremove(object = "legend")
        p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1))
        ggsave(p3, filename = paste0(figpath,"VLN ", varsall[i], ".pdf"), width = 4, height = 4)
      }
    }
  }
  # bind results
  mmres =  res_list[!is.na(res_list)] %>% dplyr::bind_rows()
  return(mmres)
}
