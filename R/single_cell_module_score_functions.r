######## weighted module score and mixed GLM for 2 group 2 timepoint repeated measures (default) comparison

WeightedCellModuleScore = function(seurat_object, module_list, threshold = 0.1, cellwise_scaling = FALSE, return_weighted = TRUE, Seurat_version = "2"){
  # init
  require(dplyr) ; require(magrittr)
  score_keep = list()
  names_vec = vector()
  library(parallel)
  #apply to each module
  print("function returns a dataframe with cell barcodes as rownames and module as column name")
  for (i in 1:length(module_list)) {
    # init
    signature = module_list[[i]]
    signaturename = names(module_list[i])
    # calc weights
    n_modules = length(module_list)

    if(Seurat_version == "2") {
      print(paste0("calculating module score for module ", i , " of ",n_modules ))
      gene_universe = intersect(signature, rownames(seurat_object@data))
      if ( (length(signature) / length(gene_universe) ) < threshold ) {
        print(paste0("< 10% of genes present in any cell for " , signaturename, " not calculating score; set threshold to 0 to override"))
        i = i + 1
      } else {
        names_vec[i] = signaturename
        # calculate the number of genes in the signature with nonzero expression in each cell
        gene_rep = apply(seurat_object@data[gene_universe, ], 2, FUN = function(x){ length(x[x>0]) / length(gene_universe) } )
        # calculate av base score to be modified by return args
        mod_avg = Matrix::colMeans(seurat_object@assays$RNA@data[gene_universe, ])
      }
    } else {
      print(paste0("calculating weighted av of module ", i , " of ", n_modules))
      gene_universe = intersect(signature, rownames(seurat_object@assays$RNA@data))
      if ((length(signature) / length(gene_universe) ) < threshold ) {
        print(paste0("< 10% of genes present in any cell for " , signaturename, " not calculating score; set threshold to 0 to override"))
        i = i + 1
      }
      names_vec[i] = signaturename
      # calculate the number of genes in the signature with nonzero expression in each cell
      gene_rep = apply(seurat_object@assays$RNA@data[gene_universe, ], 2, FUN = function(x){ length(x[x>0]) / length(gene_universe) } )
      # calculate av base score to be modified by return args
      mod_avg = Matrix::colMeans(seurat_object@assays$RNA@data[gene_universe, ])
    }
    # option to return score scaled across all cells / can change to return non weighted score
    if (isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell z scored of cells in object
      score_return = scaled_mod_weight_avg =  mod_avg * gene_rep %>% base::scale()
    }
    if (!isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell
      score_return = weighted_average = mod_avg * gene_rep
    }
    if (isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average in cell z scored of cells in object
      score_return = scaled_mod_avg =  mod_avg  %>% base::scale()
    }
    if (!isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average exprs of genes in cell
      score_return =  mod_avg
    }
    # format returned score as dataframe columns = modules, rownames = barcodes
    score_return = as.data.frame(score_return)
    names(score_return) = names_vec[i]
    score_keep[[i]] = score_return
  }
  score_df = do.call(cbind, score_keep)
  return(score_df)
  # can add this back to Seurat Object or SCE object as metadata
}


#### CALCULATE MIXED MODEL SCORE FOR 2 GROUP CONTRAST TIME DELTA AND BASELINE DELTA
# set up module data frame
# module_df = h1@meta.data %>%
#   select(barcode_check, celltype_joint, module_n)

# add metadata for lme4 model
# met = read_delim("git_ignore/full_metadata/full_sample_metadata.txt", delim = "\t")
# met = met %>% filter(CITEdata == 1) %>% filter(time_cohort == "d1" & vaccine_cohort == "H1N1")

# format metadata as factors for lme4 ordered for contrast contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
# md = h1@meta.data %>%
#   mutate(group_id = factor(adjmfc.time,  levels = c('d0 low', 'd1 low', 'd0 high', 'd1 high'))) %>%
#   mutate_if(is.character, as.factor) %>%
#   select(barcode_check, sampleid, timepoint, group_id, gender, celltype_joint)


# define contrasts for levels of group_id combined group_time factor : levels(s@meta.data$group_id) [1] "0_0" "0_1" "1_0" "1_1"
c00 = c(1,0,0,0) ; c01 = c(0,1,0,0) ; c10 = c(0,0,1,0) ; c11 = c(0,0,0,1)
contrast_2 = list("time1vs0_group2vs1" = ((c11 - c10) - (c01 - c00)),"time0_group2vs1" = (c10 - c00))
f1 = 'modulescore ~ group_id + (1|sampleid)'

# module mized model function
FitModuleMixModel = function(module_data_frame, lme4metadata, f1 = f1, contrast_list = contrast_2, plot_savepath = NULL){
  require(lme4); require(emmeans); '%ni%' = Negate('%in%')
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

# Fit mixed model
# mm_res = FitModuleMixModel(module_data_frame = module_df, lme4metadata = md, f1 = f1, contrast_list = contrast_2, plot_savepath = plot_savepath)

# module mized model function
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
