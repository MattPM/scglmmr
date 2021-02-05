#' GroupContrastGLMMsinglecell - within each cell type contrast difference in fold change between groups, baseline difference, and fold change across groups of module scores.
#'
#' @param module_data_frame data for each cell to model -- designed to be scores for modules (columns) returned by scglmmr::WeightedModuleScore
#' @param lme4metadata
#' @param f1
#' @param contrast_list
#' @param plot_savepath
#' @param celltype_column the column in metadata with the cluster / celltype designation for each cell
#'
#' @return
#' @import ggplot2
#' @import lme4
#' @importFrom stats as.formula
#' @importFrom rlang sym
#' @importFrom emmeans emmeans
#' @importFrom ggpubr ggpaired rremove
#' @importFrom egg ggarrange
#' @importFrom broom tidy
#' @importFrom dplyr select mutate bind_rows if_else filter
#' @importFrom tidyr spread everything
#' @export
#'
#' @examples
#' #### CALCULATE MIXED MODEL SCORE FOR 2 GROUP CONTRAST TIME DELTA AND BASELINE DELTA
#' # set up module data frame
#' `module_data_frame` are any data for ech cell to model -- designed to be scores for each model
#'  module_df = h1@meta.data %>%
#'   select(barcode_check, celltype_joint, module_n)
#'   # add metadata for lme4 model
#'   met = read_delim("git_ignore/full_metadata/full_sample_metadata.txt", delim = "\t")
#'  # format metadata as factors for lme4 ordered for contrast contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
# md = h1@meta.data %>%
#   mutate(group_id = factor(adjmfc.time,  levels = c('d0 low', 'd1 low', 'd0 high', 'd1 high'))) %>%
#   mutate_if(is.character, as.factor) %>%
#   select(barcode_check, sampleid, timepoint, group_id, gender, celltype_joint)
GroupContrastGLMMsinglecell = function(module_data_frame, celltype_column = 'celltype', metadata,
                                        fixed_effects = NULL, lmer_formula = NULL, plotdatqc = TRUE, figpath){


  # specify custom contrasts difference in treatment between groups, treatment effect across groups, per-treatment difference between groups.
  c00 = c(1,0,0,0)
  c01 = c(0,1,0,0)
  c10 = c(0,0,1,0)
  c11 = c(0,0,0,1)
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
    f1 = paste("modulescore ~ 0 +", fixef, '+ group_id + (1|subjectid)')
  } else {
    f1 = paste("modulescore ~ 0 +", 'group_id + (1|subjectid)')
  }
  f_store = f1
  f1 = stats::as.formula(f1)

  # checks on minimal required metadata and formatting of group levels
  print('model specified:'); print(f1)
  stopifnot(is.factor(metadata$group_id))
  cat("check levels, testing: \n",
      "1: difference in treatmant effect between groups: ",
      levels(metadata$group_id)[3:4], " fold change  vs: ",
      levels(metadata$group_id)[1:2], "fold change  \n",
      "2: pre treatment baseline difference between groups: ",
      levels(metadata$group_id)[c(1,3)], " fold change \n",
      "3: treatment effect across both groups combined \n",
      "a prori contrasts across `ggroup_id` variable: "
  )
  print(contrast_list)

  # print data format note
  # print("module_data_frame columns: 1 = cell barcode 2 = celltype 3-nmodules = name of each module")
  #print("celltype column of lme4metadata and module_data_frame mmust be 'celltype_joint' ")

  # init storage
  res_celltype = res_list = list()

  # define model data
  module_names = names(module_data_frame)
  cts = as.character(unique(metadata[[celltype_column]]))

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))

    # subset metadat to dataeeded for modeling

    metsub = metadata[metadata[[celltype_column]] == cts[u], ]
    meta_fit = colnames(metsub)
    stopifnot(c( 'group_id', 'subjectid') %in% meta_fit)

    # subset module dataframe to cells in celltype u from metadata above
    mod_df = module_data_frame[rownames(metsub), ]

    # merge module data and model metadata
    dat = cbind(metsub, mod_df)

    # indexing over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]; print(module)
      dat_fit = dat[ ,c(module, meta_fit)]
      colnames(dat_fit)[1] = "modulescore"

      # fit models
      gvar = rlang::sym('group_id')
      m1 = tryCatch(lme4::lmer(formula = f1, data = dat_fit), error = function(e)  return(NA))
      emm1 = tryCatch(emmeans::emmeans(object = m1, specs = ~ {{ group_id }}, data = dat_fit, lmer.df = "asymptotic"), error = function(e) return(NA))

      # error checking
      if( suppressWarnings(is.na(m1))) {
        print(paste0(" could not fit  model for ", module_names[i]))
        res_list[[i]] = NA
      }
      if( suppressWarnings(is.na(emm1))) {
        print(paste0(module_names[i],  " check vcov matrix could not calculate marginal means  "))
        res_list[[i]] = NA
      } else{

        # format model results
        marginal_mean_df = broom::tidy(emm1)[ ,1:2] %>% tidyr::spread('group_id','estimate')
        names(marginal_mean_df) = paste0(names(marginal_mean_df), "_marginal_mean")

        # apply custom contrast
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
        contrastdf = contrastdf %>%
          dplyr::mutate(singular_fit = dplyr::if_else(is.null(singular_fit), true = 0, false = 1))


        # add a Wilcoxon test of baseline difference between groups
        baseline_data = dat_fit %>% filter(group_id %in% levels(dat_fit$group_id)[c(1,3)])
        baseline_data$group_id = droplevels(baseline_data$group_id)
        wt = wilcox.test(modulescore ~ group_id, data = baseline_data) %>%
          broom::tidy() %>%
          rename(p.value.t0.wilcox = p.value)

        # merge all results and  format into a single line of a list to be r bound.
        contrastdf = cbind(contrastdf, wt)

        if(isTRUE(plotdatqc)){
          # reformat minimal marginal means
          p1 = plot(emm1) +
            coord_flip() +
            theme_classic() +
            theme(plot.title = element_blank()) +
            xlab("model estimated marginal means") +
            theme(axis.title.y = element_text(size = 7)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), axis.title.x = element_blank()) +
            scale_x_continuous(position = "left") +
            theme(axis.text.y = element_text(size = 4)) +
            theme(axis.line = element_line(colour = 'black', size = 0.2)) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


          ## add baseline and FC p value to plot
          t0p = formatC(contrastdf$p.value.t0.wilcox, format = "e", digits = 4)
          t1p = formatC(contrastdf$p.valuetime1vs0_group2vs1, format = "e", digits = 4)
          titleplot = paste0(cts[u]," marginal means ",module)

          p2 = ggplot(dat_fit, aes(x = group_id, y = modulescore)) +
            geom_violin(aes(fill = group_id), show.legend = FALSE) +
            # geom_boxplot(width=0.14, outlier.shape = NA) +
            stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", color = "black", size = 0.3) +
            scale_fill_manual(values = c("#1E90FFD9", "#191970D9", "#FF0000D9", "#CD0000D9")) +
            ggtitle(paste0(f_store, "\n", cts[u], " ", module, "\n", "baseline wilcox p = ", t0p, "\n",
                           "lmer 24h Fold Change delta p = ", t1p)) +
            theme_bw(base_size = 11) +
            theme(axis.title.y =  element_text(face = "bold", size = 7)) +
            theme(plot.title = element_text(size = 7, face = "bold")) +
            ylab(paste0(module, "  score  " )) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7, face = "bold")) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

          # combine marginal means and data view
          p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1))
          ggsave(p3, filename = paste0(figpath,"VLN ", titleplot, ".pdf"), width = 3.5, height = 3.5)
        }
        # store result vector in result list i for celltype u
        res_list[[i]] =
          cbind(contrastdf, marginal_mean_df) %>%
          dplyr::mutate(module = module) %>%
          dplyr::mutate(formula = f_store) %>%
          dplyr::select(module, tidyr::everything())
      }
    }
    res_list = res_list[!is.na(res_list)]
    resdf = do.call(rbind, res_list)

    # store bound results for each module for celltype u in list element u
    res_celltype[[u]] = cbind(celltype = cts[u], resdf)
  }
  # bind results from each cell type into single dataframe
  resdf_full = do.call(rbind, res_celltype)
  return(resdf_full)
}


