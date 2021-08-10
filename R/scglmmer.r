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
#f1 = 'modulescore ~ group_id + (1|sampleid)'

# module mized model function
FitGLM = function(gene_data_frame, lme4metadata, f1 = f1, contrast_list = contrast_2, plot_savepath = NULL){
  #require(lme4); require(emmeans);
  '%ni%' = Negate('%in%')
  print(paste0("model currently specified = poisson GLMM of ", f1, " change argument to f1 to add factor variables in lme4metadata argument"))
  # minimal required metadata
  stopifnot(is.factor(lme4metadata$group_id)) ; stopifnot(is.factor(lme4metadata$sampleid))
  # print data format note
  print("group_id level1 must be group 1 baseline (intercept) e.g. md %>% mutate(group_id = factor(group_timepoint,  levels = c('0_0', '0_1', '1_0', '1_1')))")
  print("gene_data_frame columns: 1 = cell barcode 2 = celltype 3-nmodules = name of each Gene ")
  print("celltype column of lme4metadata and gene data gene_data_frame must be 'celltype_joint' ")

  # Test
  # gene_data_frame = module_df
  # lme4metadata = md
  # f1 = f1
  # contrast_list = contrast_2
  # plot_savepath = plot_savepath
  # cts = unique(gene_data_frame$celltype_joint)[c(1,2)]

  # init storage
  res_celltype = res_list = list()

  # define module names
  module_names = names(gene_data_frame)
  module_names = module_names[-c(1:2)]

  # define celltypes
  cts = unique(gene_data_frame$celltype_joint) %>% as.character()

  # show what is being tested
  print("testing the following genes in these clusters and comparing ") ; print(names(contrast_list))
  print(contrast_list);  print(cts); print(module_names)

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())

    # celltype indexed
    print(paste0("fitting models for ", cts[u]))
    mod_df = gene_data_frame %>% filter(celltype_joint == cts[u]) %>% select(-celltype_joint)
    met = lme4metadata %>% filter(celltype_joint %in% cts[u])
    # definee columns of lme4 metadata for use in next section; combine meta and module score
    meta_fit = colnames(lme4metadata)
    dat = full_join(x = met, y = mod_df, by = "barcode_check")
    # indexed over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]; print(module)
      dat_fit = dat %>% select(modulescore = module, meta_fit)

      # fit model and calculate margina, meads
      m1 = tryCatch(lme4::glmer(f1, data = dat_fit, family = poisson()), error = function(e)  return(NA))
      emm1 = tryCatch(emmeans::emmeans(object = m1, specs = ~ group_id, data = dat_fit), error = function(e) return(NA))

      # optional plot save of means
      if( !is.null(plot_savepath) & !is.na(emm1) ){
        titleplot = paste0(cts[u]," marginal means ",module)
        tryCatch(plot(emm1, comparisons = TRUE) +
                   ggtitle(titleplot) + theme_bw(base_size = 18) +
                   theme(plot.title = element_text(size = 8)) +
                   coord_flip() +
                   ggsave(filename = paste0(plot_savepath, titleplot, ".png"), width = 4.5, height = 5),
                 error = function(e) e )
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
# mm_res = FitModuleMixModel(gene_data_frame = module_df, lme4metadata = md, f1 = f1, contrast_list = contrast_2, plot_savepath = plot_savepath)


