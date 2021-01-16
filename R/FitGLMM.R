# poisson random intercept model on counts
# specify model formulae 
f1 = 'gene ~ offset(log(nUMI)) + timepoint + (1 | sampleid)'
f2 = 'gene ~ offset(log(nUMI)) +  (1 | sampleid)'

FitGLMM = function(gene_data, lme4metadata, reduced_model = f2, full_model = f1, celltype_genes_test = NULL){ 
  require(lme4); require(lmtest); '%ni%' = Negate('%in%')
  print(paste0("full model currently specified = poisson GLMM of ", f1, " vs reduced model ", f2))
  #print(paste0("only fit random interept model. Random slope models not supported (effect size will be incorrect) "))
  
  # # Test function 
  #module_data_frame = gene_data
  #lme4metadata = md
  #reduced_model = f2
  #full_model = f1
  #celltype_genes_test = genes_test
  #f1 = 'gene ~ offset(log(nUMI)) + timepoint + (1 | sampleid)'
  #f2 = 'gene ~ offset(log(nUMI)) +  (1 | sampleid)'
  
  ## Run checks for minimal required metadata
  if (c("celltype") %ni% colnames(lme4metadata)) {
    stop("in the argument lme4metadata, please rename your cluster column 'celltype' this is minimal required metadata for the model ")
  }
  if (c("barcode") %ni% colnames(lme4metadata)) {
    stop("in the argument lme4metadata, please rename your cell barcode column 'barcode' this is required for internal function subsetting ")
  }
  if (c("nUMI") %ni% colnames(lme4metadata)) {
    stop("in the argument lme4metadata, please rename your number of umi 'nUMI' this is minimal required metadata for the model ")
  }
  
  # Run checks for offset term that should be specified in Poisson model 
  if ( isFALSE(grepl("nUMI" , f1, fixed = TRUE))) {
    stop("to adjust for exposure bias a NUMI offset should be used in the model ")
  }
  if ( isFALSE(grepl("nUMI" , f2, fixed = TRUE))) {
    stop("to adjust for exposure bias a nUMI offset should be used in the model formula, e.g. f1 =  ")
  }
  
  # Run checks for minimal required structure of f1 and f2 lme4 formulae
  if ( isFALSE(grepl("gene" , f1, fixed = TRUE))) {
    stop(" 'gene' not found in f1, the formula should be structured expr ~ offset(log(nUMI)) + your_fixed_effects + (1 | your_varying_effects ) ")
  }
  if ( isFALSE(grepl("gene" , f2, fixed = TRUE))) {
    stop(" 'gene' not found in f2, the formula should be structured expr ~ offset(log(nUMI)) + your_fixed_effects + (1 | your_varying_effects ) ")
  }
  
  # init storage 
  res_celltype = res_list = list() 
  
  # define gene names to test 
  gene_names = colnames(gene_data)
  #gene_names = gene_names[-c(1:2)]
  
  # define celltypes
  cts = unique(md$celltype) %>% as.character()
  
  # make sure celltypes ordered by subset of genes test list by cell type   
  if(!is.null(celltype_genes_test)){
    stopifnot(length(cts) == length(celltype_genes_test))
    cts = names(celltype_genes_test)
  } 
  

  # 1. Indexed over cell type: 
  for (u in 1:length(cts)) {
    suppressMessages(gc()) # garbage collect bc of as matrix call 
    
    # celltype indexed 
    met_ = lme4metadata %>% filter(celltype %in% cts[u])
    df_ = gene_data[met_$barcode, ]
    print(paste0("fitting models for ", unique(met_$celltype)))
    
    # define columns of lme4 metadata for use in next section; combine meta and module score 
    #meta_fit = colnames(lme4metadata)
    dat = cbind(met_, as.data.frame(as.matrix(df_))) 
    
    #dat = full_join(x = met_ , y = df_ , by = "barcode") 
    
    if( ! is.null(celltype_genes_test)){
      dat = dat %>% select(colnames(lme4metadata), celltype_genes_test[[u]])
      gene_names = celltype_genes_test[[u]]
      print(paste0("fitting mixed model on ",  length(gene_names), " genes" ))
    } 
    
    # 2. Indexed over gene
    for (i in 1:length(gene_names)) {
      # data for gene i 
      gene_ = gene_names[i]; print(gene_) 
      dat_fit = dat %>% select(gene = gene_, colnames(lme4metadata))
     # dat_fit = as.matrix(dat_fit)
      # fit models
      m1 = tryCatch(lme4::glmer(f1, data = dat_fit, family = poisson(link = "log")), error = function(e) return(NA))
      m2 = tryCatch(lme4::glmer(f2, data = dat_fit, family = poisson(link = "log")), error = function(e) return(NA))
      
      #Error handling 
      if( suppressWarnings(is.na(m1))) { 
        print(paste0(gene_names[i],  " could not fit  model for ", as.character(cts[u])))
        res_list[[i]] = NA 
      }
      if(suppressWarnings(is.na(m2))) { 
        print(paste0(gene_names[i],  " could not fit model for  ", as.character(cts[u])))
        res_list[[i]] = NA 
      } else{  
        
        # likelihood ratio test  ; fixed_effect = coef(m1)[[1]][1 ,2] all coefficients are the same fixed slope model
        lrtest = lrtest(m1, m2)
        pval = lrtest$`Pr(>Chisq)`[2]
        fixed_effect = coef(m1)[[1]][1 ,2]
        # flag singular fit 
        singular_fit = m2@optinfo$conv$lme4$messages
        # return results 
        result = bind_cols(p.value = pval, fixed_effect = fixed_effect, singularity = singular_fit) %>% 
          mutate(singular_fit = if_else(is.null(singular_fit), true = 0, false = 1))
        
        # merge into result list 
        res_list[[i]] = result %>%
          mutate(gene = gene_) %>% 
          mutate(fullmodel = f1) %>% 
          mutate(reducedmodel = f2) %>% 
          select(gene, everything())
      }
    }
    # remove all module / celltypes without model fit; bind to df; add to cellype list. 
    #res_list = res_list[!is.na(res_list)]
    resdf = bind_rows(res_list)
    res_list = list()
    # indexed over celltype store celltype 1 values in resdfct[[u]]  
    res_celltype[[u]] = resdf %>% mutate(celltype = cts[u]) %>% 
      mutate(padj = p.adjust(p = p.value, method = "BH")) %>% 
      select(celltype, gene, p.value, padj, everything())
  }
  resdf_full = bind_rows(res_celltype)
  return(resdf_full)
}
