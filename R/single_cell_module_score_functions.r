
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
  c00 = c(1,0,0,0) ; c01 = c(0,1,0,0) ; c10 = c(0,0,1,0) ; c11 = c(0,0,0,1)
  contrast_list = list(
    "time1vs0_group2vs1" = (c11 - c10) - (c01 - c00),
    "time1vs0" = (c11 + c01) / 2 - (c10 + c00) / 2,
    "time0_group2vs1" = c10 - c00
  )

  ### testfun
  lmer_formula = NULL
  module_data_frame = mdf
  metadata = met
  celltype_column = "celltype_joint"
  fixed_effects = 'batch'


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

  # define module names
  module_names = names(module_data_frame)
  cts = as.character(unique(metadata[[celltype_column]]))

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))

    # subset metadat to dataeeded for modeling
    metsub = met[met[[celltype_column]] == cts[u], c(fixed_effects, 'group_id', 'subjectid')]
    meta_fit = colnames(metsub)

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
        print(paste0(" could not fit  model for ", varsall[i]))
        res_list[[i]] = NA
      }
      if(is.na(emm1)){
        print(paste0(varsall[i],  " check vcov matrix could not calculate marginal means  "))
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
        res_list[[i]] =
          cbind(contrastdf, marginal_mean_df) %>%
          dplyr::mutate(module = module) %>%
          dplyr::mutate(formula = f_store) %>%
          dplyr::select(module, tidyr::everything())

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
            scale_fill_manual(values = c("dodgerblue", "midnightblue", "red", "firebrick")) +
            ggtitle(paste0(f_store, "\n", cts[u], " ", module, "\n", "baseline wilcox p = ", t0p, "\n",
                           "lmer 24h Fold Change delta p = ", t1p)) +
            theme_bw(base_size = 11) +
            theme(axis.title.y =  element_text(face = "bold", size = 7)) +
            theme(plot.title = element_text(size = 7, face = "bold")) +
            ylab(paste0(module, "  score  " )) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

          # combine marginal means and data view
          p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1))
          ggsave(p3, filename = paste0(figpath,"VLN ", titleplot, ".pdf"), width = 3.5, height = 3.5)
        }
      }
    }
    # remove all module / celltypes without model fit; bind to df; add to cellype list.
    res_list = res_list[!is.na(res_list)]
    resdf = do.call(rbind, res_list)
    resdf = cbind(celltype = cts[u], resdf)
  }
  resdf_full = do.call(rbind, res_celltype)
  return(resdf_full)
}




# Fit mixed model
# mm_res = FitModuleMixModel(module_data_frame = module_df, lme4metadata = md, f1 = f1, contrast_list = contrast_2, plot_savepath = plot_savepath)
# module mized model function
#' WilcoxWithinCluster
#'
#' @param module_data_frame
#' @param lme4metadata
#' @param plot_savepath
#'
#' @return
#' @export
#'
#' @examples
WilcoxWithinCluster = function(module_data_frame, lme4metadata, plot_savepath = NULL){
  # require(lme4); require(emmeans);
  '%ni%' = Negate('%in%')
  # emm_options(lmerTest.limit = Inf)
  # emm_options(pbkrtest.limit = Inf)

  # minimal required metadata
  #stopifnot(is.factor(lme4metadata$group_id)) ;
  #stopifnot(is.factor(lme4metadata$sampleid))
  # print data format note
  # print("group_id level1 must be group 1 baseline (intercept) e.g. md %>% mutate(group_id = factor(group_timepoint,  levels = c('0_0', '0_1', '1_0', '1_1')))")
  print("module_data_frame columns: 1 = cell barcode 2 = celltype 3-nmodules = name of each module")
  print("celltype column of lme4metadata and module_data_frame mmust be 'celltype_joint' ")

  # Test
  #  module_data_frame = m2
  #lme4metadata = md_colitis_vs_nocolitis
  # f1 = f1
  # contrast_list = contrast_2
  # plot_savepath = plot_savepath
  # cts = unique(module_data_frame$celltype_joint)[c(1,2)]

  # init storage
  res_celltype = res_list = list()

  # define module names
  module_names = names(module_data_frame)
  module_names = module_names[-c(1:2)]

  # define celltypes
  cts = unique(module_data_frame$celltype_joint) %>% as.character()

  # show what is being tested
  print("testing the following modules in these clusters and comparing ")
  #print(contrast_list);
  print(cts); print(module_names)

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))
    mod_df = module_data_frame %>% filter(celltype_joint == cts[u]) %>% select(-celltype_joint)
    met = lme4metadata %>% filter(celltype_joint %in% cts[u])
    # definee columns of lme4 metadata for use in next section; combine meta and module score
    meta_fit = colnames(lme4metadata)
    dat = full_join(x = met, y = mod_df, by = "barcode_check")
    # indexed over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]; print(module)
      dat_fit = dat %>% select(modulescore = module, meta_fit)
      # baseline simple DE test with paired wilcox test
      baseline_data = dat_fit

      # %>% filter(group_id %in% levels(dat_fit$group_id[c(1,3)]))
      wt = wilcox.test(modulescore ~ group_id, data = baseline_data
                       #%>% filter(group_id %in% levels(dat_fit$group_id)[c(1,3)])
      ) %>%
        broom::tidy() %>%
        rename(p.value.t0.wilcox = p.value)

      # contrast df
      # contrastdf = cbind(contrastdf, wt)
      contrastdf = wt
      t0p = formatC(contrastdf$p.value.t0.wilcox, format = "e", digits = 4)

      # plot marginal and violin
      if( !is.null(plot_savepath) ){

        # reformat minimal marginal means

        titleplot = paste0(cts[u],module)
        p3 = ggplot(baseline_data, aes(x = group_id, y = modulescore, fill = group_id)) +
          geom_jitter(shape = 16, alpha = 0.5, size = 0.5, show.legend = FALSE) +
          geom_violin(show.legend = FALSE) +
          scale_fill_manual(values = c("darkorange2", "dodgerblue3"))+
          theme_bw() +
          theme(axis.title.y = element_text(size = 13)) +
          theme(axis.text.x = element_text(size = 13), axis.title.x = element_blank()) +
          theme(axis.text.y = element_text(size = 13)) +
          theme(axis.line = element_line(colour = 'black', size = 0.2)) +
          ggtitle(paste0(cts[u], " ", module, "\n", "baseline wilcox p = ", t0p))
        ggsave(p3, filename = paste0(plot_savepath,"VLN ", titleplot, ".pdf"), width = 3.3, height = 3.3)
        ggsave(p3, filename = paste0(plot_savepath,"VLN ", titleplot, ".png"), width = 3.3, height = 3.3)
        #ggsave(p3, filename = paste0(testpath, "test2.pdf"), width = 3, height = 3)

      }

      # merge into result list
      res_list[[i]] =
        cbind(contrastdf) %>%
        mutate(modulename = module) %>%
        select(modulename, everything())
    }

    # remove all module / celltypes without model fit; bind to df; add to cellype list.
    res_list = res_list[!is.na(res_list)]
    resdf = bind_rows(res_list)
    # indexed over celltype store celltype 1 values in resdfct[[u]]
    res_celltype[[u]] = resdf %>% mutate(celltype = cts[u]) %>% select(celltype, everything())
  }
  resdf_full = bind_rows(res_celltype)
  return(resdf_full)
}



#### Depreicted functions
FitModuleMixModel = function(module_data_frame, lme4metadata, f1 = f1, contrast_list = contrast_2, plot_savepath = NULL){
  .Deprecated("GroupContrastGLMMsinglecell")
  #require(lme4); require(emmeans); '%ni%' = Negate('%in%')


  # emm_options(lmerTest.limit = Inf)
  # emm_options(pbkrtest.limit = Inf)
  print(paste0("model currently specified  = ", f1))
  # minimal required metadata
  stopifnot(is.factor(lme4metadata$group_id)) ; stopifnot(is.factor(lme4metadata$sampleid))
  # print data format note
  print("group_id level1 must be group 1 baseline (intercept) e.g. md %>% mutate(group_id = factor(group_timepoint,  levels = c('0_0', '0_1', '1_0', '1_1')))")
  print("module_data_frame columns: 1 = cell barcode 2 = celltype 3-nmodules = name of each module")
  print("celltype column of lme4metadata and module_data_frame mmust be 'celltype_joint' ")

  # Test
  # module_data_frame = module_df
  # lme4metadata = md
  # f1 = f1
  # contrast_list = contrast_2
  # plot_savepath = plot_savepath
  # cts = unique(module_data_frame$celltype_joint)[c(1,2)]

  # init storage
  res_celltype = res_list = list()

  # define module names
  module_names = names(module_data_frame)
  module_names = module_names[-c(1:2)]

  # define celltypes
  cts = unique(module_data_frame$celltype_joint) %>% as.character()

  # show what is being tested
  print("testing the following modules in these clusters and comparing ") ; print(names(contrast_list))
  print(contrast_list);  print(cts); print(module_names)

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))
    mod_df = module_data_frame %>% filter(celltype_joint == cts[u]) %>% select(-celltype_joint)
    met = lme4metadata %>% filter(celltype_joint %in% cts[u])
    # definee columns of lme4 metadata for use in next section; combine meta and module score
    meta_fit = colnames(lme4metadata)
    dat = full_join(x = met, y = mod_df, by = "barcode_check")
    # indexed over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]; print(module)
      dat_fit = dat %>% select(modulescore = module, meta_fit)

      #dat = cbind(modulescore = module_dataframe[ ,module], lme4metadata) %>% as.data.frame()
      m1 = tryCatch(lmer(formula = f1, data = dat_fit), error = function(e)  return(NA))
      emm1 = tryCatch(emmeans(object = m1, specs = ~ group_id, data = dat_fit, lmer.df = "asymptotic"), error = function(e) return(NA))

      # optional plot save of means
      if( !is.null(plot_savepath) & !is.na(emm1) ){
        titleplot = paste0(cts[u]," marginal means ",module)
        tryCatch(
          plot(emm1) +
            ggtitle(paste0(titleplot, "\n", f1)) +
            theme_bw(base_size = 11) +
            theme(plot.title = element_text(size = 8, face = "bold")) +
            xlab("Model estimated marginal means adjusted for covariates") +
            coord_flip() +
            ggsave(filename = paste0(plot_savepath, titleplot, ".pdf"), width = 4.5, height = 5),
          error = function(e) e
        )

      }
      #Error handling
      if( suppressWarnings(is.na(m1))) {
        print(paste0(module_names[i],  " could not fit  model for ", as.character(cts[u])))
        res_list[[i]] = NA
      }
      if(is.na(emm1)){
        print(paste0(module_names[i],  " check cov matrix could not calculate marginal means for  ", as.character(cts[u])))
        res_list[[i]] = NA
      } else{
        marginal_mean_df = broom::tidy(emm1) %>% select(group_id, estimate) %>%  spread(group_id, estimate)
        names(marginal_mean_df) = paste0(names(marginal_mean_df), "_marginal_mean")
        # specify contrasts and fit contrasts
        contrast_fit = contrast(emm1, method = contrast_list)
        # format dataframe from model results
        contrastdf = broom::tidy(contrast_fit)
        contrast_names = names(contrast_list)
        contrast_colnames = sapply(contrast_names, function(x) paste0(colnames(contrastdf), x)) %>% as.data.frame() %>% unlist(use.names = F) %>% as.character()
        contrastdf = cbind(contrastdf[1, ], contrastdf[2, ]); names(contrastdf) = contrast_colnames

        # flag singular fit
        singular_fit = m1@optinfo$conv$lme4$messages
        contrastdf = contrastdf %>% mutate(singular_fit = if_else(is.null(singular_fit), true = 0, false = 1))

        # baseline simple DE test with paired wilcox test
        baseline_data = dat_fit %>% filter(group_id %in% levels(dat_fit$group_id[c(1,3)]))
        wt = wilcox.test(modulescore ~ group_id, data = baseline_data %>% filter(group_id %in% levels(dat_fit$group_id)[c(1,3)]) ) %>%
          broom::tidy() %>%
          rename(p.value.t0.wilcox = p.value)

        # merge
        contrastdf = cbind(contrastdf, wt)


        # plot marginal and violin
        if( !is.null(plot_savepath) & !is.na(emm1) ){

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
            scale_fill_manual(values = c("dodgerblue", "midnightblue", "red", "firebrick")) +
            ggtitle(paste0(f1, "\n", cts[u], " ", module, "\n", "baseline wilcox p = ", t0p, "\n",
                           "lmer 24h Fold Change delta p = ", t1p)) +
            theme_bw(base_size = 11) +
            theme(axis.title.y =  element_text(face = "bold", size = 7)) +
            theme(plot.title = element_text(size = 7, face = "bold")) +
            ylab(paste0(module, "  score  " )) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

          # combine marginal means and data view
          p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1))
          ggsave(p3, filename = paste0(plot_savepath,"VLN ", titleplot, ".pdf"), width = 3, height = 3)
          #ggsave(p3, filename = paste0(testpath, "test2.pdf"), width = 3, height = 3)

          # save qqlpot
          # plot model quality control
          pdf(paste0(plot_savepath,"qc ", titleplot,".pdf"), width = 9, height = 9)
          p4 = lattice::qqmath(m1)
          p5 = lattice::levelplot(as.matrix(vcov.merMod(m1)),
                                  main = "variance covariance matrix",
                                  scales=list(y=list(cex=.5),x=list(rot=90, cex=.5)),
                                  col.regions = viridis::plasma(100) )
          p6 = lattice::dotplot(ranef(m1))
          p7 = qqmath(ranef(m1))
          gridExtra::grid.arrange(p4, p5, p6$sampleid, p7$sampleid)
          dev.off()

        }

        # merge into result list
        res_list[[i]] =
          cbind(contrastdf, marginal_mean_df) %>%
          mutate(modulename = module) %>%
          mutate(formula = f1) %>%
          select(modulename, everything())
      }
    }
    # remove all module / celltypes without model fit; bind to df; add to cellype list.
    res_list = res_list[!is.na(res_list)]
    resdf = bind_rows(res_list)

    # indexed over celltype store celltype 1 values in resdfct[[u]]
    res_celltype[[u]] = resdf %>% mutate(celltype = cts[u]) %>% select(celltype, everything())
  }
  resdf_full = bind_rows(res_celltype)
  return(resdf_full)
}
