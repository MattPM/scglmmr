#' GroupContrastGLMMtidy - fit a random intercept model of the difference in fold change between groups and fold change across both groups and baseline controlling for any number of covariates.
#'
#' @param tidydf tidy data must contain column 'subjectid' and any names covariates (specified with )
#' @param response_variable the column containing the variables to be modeled in tidy (long) format -- e.g. 'gene' with the gene's values in the column in response_value
#' @param response_value the value for the measurement response_variable -- should probably be in log space.
#' @param lmer_formula lmer model formula, if not specified, creates a donor random intercept model of 0 + group_id + (1|subjectid)
#' @param figpath path to save figures e.g.  file.path()
#' @param fixed_effects the fixed effects to specify in this case as covariates e.g. c('gender', 'ethnicity'), these are automatically combined into formula with your_fixed_effects + 0 + group_id + (1|subjectid)
#' @return
#' @import ggplot2
#' @importFrom lme4 lmer vcov.merMod
#' @importFrom emmeans emmeans
#' @importFrom lattice qqmath levelplot dotplot
#' @importFrom gridExtra grid.arrange
#' @importFrom viridis plasma
#' @importFrom ggpubr ggpaired rremove
#' @importFrom egg ggarrange
#' @importFrom broom tidy
#' @importFrom dplyr select mutate bind_rows if_else
#' @importFrom tidyr spread everything
#' @export
#'
#' @examples
#' \donttest{
#'   tidydf = d2 %>%
#'   dplyr::rename('subjectid'= subject) %>%
#'   mutate(group_id = factor(group_id, levels = c(0_0, 0_1, 1_0, 1_1)))
#'
#'   response_variable = "logprot"
#'   fixed_effects = c("Age", "Race", "Gender")
#'   mmres = GroupContrastGLMMtidy(tidydf = tidydf, response_variable = 'protein', response_value = 'log10protvalue', fixed_effects = c('gender', 'ethnicity'), lmer_formula = NULL, subject_colum 'subjectid', figpath = '~' )
#'   }
#'
GroupContrastGLMMtidy = function(tidydf, response_variable, response_value, fixed_effects = NULL, lmer_formula = NULL, figpath){

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

  # init store
  res_list = list()
  varsall = unique(tidydf[[response_variable]])
  # fit model for each protein
  for (i in 1:length(varsall)) {

    # fit model separately for each protein
    d3 = tidydf[tidydf[[response_variable]] == varsall[i], ]
    lme4::lmer(formula = f1, data = d3, REML = TRUE)
    m1 = tryCatch(lme4::lmer(formula = f1, data = d3, REML = TRUE), error = function(e)  return(NA))
    emm1 = tryCatch(emmeans::emmeans(m1, specs = ~group_id), error = function(e)  return(NA))

    # error checking
    if( suppressWarnings(is.na(m1))) {
      print(paste0(" could not fit  model for ", varsall[i]))
      res_list[[i]] = NA
    }
    if(is.na(emm1)){
      print(paste0(varsall[i],  " check cov matrix could not calculate marginal means  "))
      res_list[[i]] = NA
    } else{

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
        ggpubr::ggpaired(d3,x = 'group_id', y = response_value, fill = 'group_id') +
        scale_fill_manual(values = c("dodgerblue", "midnightblue", "red", "firebrick")) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90)) +
        ggtitle(varsall[i]) +
        ylab("modeled value") +
        xlab("group_id") +
        ggpubr::rremove(object = "legend")
      p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1))
      ggsave(p3, filename = paste0(figpath,"VLN ", proteins[i], ".pdf"), width = 4, height = 4)

      # format model results
      marginal_mean_df =
        broom::tidy(emm1) %>%
        dplyr::select(group_id, estimate) %>%
        tidyr::spread(group_id, estimate)
      names(marginal_mean_df) = paste0(names(marginal_mean_df), "_marginal_mean")

      contrast_fit = emmeans::contrast(emm1, method = contrast_list)
      contrastdf = broom::tidy(contrast_fit)
      contrast_names = names(contrast_list)
      contrast_colnames = sapply(contrast_names, function(x) paste0(colnames(contrastdf), x)) %>%
        as.data.frame() %>%
        unlist(use.names = F) %>%
        as.character()


      contrastdf = cbind(contrastdf[1, ], contrastdf[2, ], contrastdf[3, ])
      names(contrastdf) = contrast_colnames

      # flag singular fits
      singular_fit = m1@optinfo$conv$lme4$messages
      contrastdf = contrastdf %>% dplyr::mutate(singular_fit = dplyr::if_else(is.null(singular_fit), true = 0, false = 1))

      # merge models
      res_list[[i]] =
        cbind(contrastdf, marginal_mean_df) %>%
        dplyr::mutate(variable = varsall[i]) %>%
        dplyr::mutate(formula = f1) %>%
        dplyr::select(variable, tidyr::everything())
    }
  }
  # bind results
  mmres =  res_list[!is.na(res_list)] %>% dplyr::bind_rows()
  return(mmres)
}
