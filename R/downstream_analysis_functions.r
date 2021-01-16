#### De functions not dependent on Seurat versions. Downstream enrichment etc. 
# from 

select = dplyr::select
'%ni%' = Negate('%in%')

# https://github.com/MattPM/de_workflow/blob/master/downstream_analysis_functions.r
# on date : April 30, 2020 
# Extract a table of effect sizes and p values from the model fit; works on one model and contrast -- 
GetContrastResults = function(limma.fit.object.list, sort.by = "p", coefficient.number, contrast.name, celltypes.vector = celltypes){
  print("if using lme4 with random effects (limma-voom-lme4 implemented in dream via variancePartition) do not use eBayes, 
        instead use *GetContrastResultsRaw* directly on the fit object returned by dream https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  require(limma)
  test = list()
  for (i in 1:length(limma.fit.object.list)) {
  
    test[[i]] = 
      topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf,sort.by = sort.by) %>% 
      rownames_to_column("gene") %>% 
      mutate(contrast = rep(contrast.name)) %>% 
      mutate(pvalue0.01 = if_else(P.Value < 0.01, true = "1", false = "0")) %>% 
      mutate(celltype = celltypes.vector[i])
  }
  names(test) = celltypes.vector
  return(test)
}

GetContrastResultsRaw =  function(limma.fit.object.list, sort.by = "p", coefficient.number, contrast.name, celltypes.vector = celltypes){
  print("returning lme4 random effects raw and FDR corrected p values with dream estimated avg exp and logFC if no random effects in model run eBayes then Get ContrastResults")
  ## pt 1 return ONLY gene, logFC, AveExpr, contrast, celltype from eBayes call to format data and add raw p values from lme4, correct with BH. 
  # p value of model fit 
  pvals = lapply(limma.fit.object.list, function(x){ x$pValue[ , coefficient.number] %>% as.data.frame() %>% rownames_to_column("gene") %>% rename(P.Value = '.') })
  
  # parameters to run ordinary t statistic 
  coef = lapply(limma.fit.object.list, function(x){ x$coefficients[ ,coefficient.number] })
  stdev = lapply(limma.fit.object.list, function(x){ x$stdev.unscaled[ ,coefficient.number] })
  sigma_ = lapply(limma.fit.object.list, function(x){ x$sigma })
  # get results 
  ebf = lapply(limma.fit.object.list, eBayes) ## note eBayes t statistics are NOT used, only for formatting output, next section adds raw p values and logFC returned by lme4 
  contrast_result = GetContrastResults(limma.fit.object.list = ebf, coefficient.number = coefficient.number, sort.by = "p", contrast.name = contrast.name, celltypes.vector = celltypes.vector) 
  contrast_result = lapply(contrast_result, function(x){ x %>% select(gene, logFC, AveExpr, contrast, celltype) })
  result = t_statistic = list()
  for (i in 1:length(limma.fit.object.list)) {
    # calculate ordinary t statistic
    # Ordinary t-statistic (from limma eBayes documentation) ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
    t_statistic[[i]] = coef[[i]] / stdev[[i]] / sigma_[[i]] %>% as.data.frame()
    t_statistic[[i]] =  t_statistic[[i]] %>% rownames_to_column("gene") 
    colnames(t_statistic[[i]]) = c("gene", "t")
    
    result[[i]] = 
      contrast_result[[i]] %>% 
      mutate(P.Value = plyr::mapvalues(x = gene, from = pvals[[i]]$gene, to = round(pvals[[i]]$P.Value, digits = 9))) %>% 
      mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH")) 
    # add t statistic 
    result[[i]] = full_join(result[[i]], t_statistic[[i]], by = "gene")
    
  }
  names(result) = celltypes.vector
  return(result)
}
#Usage 
# d1hvl_res = get.contrast.results(limma.fit.object.list = ebf, coefficient.number = 1, contrast.name = "day_1_highvslow",celltypes.vector = celltypes)

# get a gene matric for plotting a statistic (logFC) from pseudobulk model results. 
GetGeneMatrix = function(result.list, gene_subset, stat_for_matrix = "logFC", pvalfilter, logfcfilter){ 
  mtx = lapply(result.list, function(x){ 
    x = x %>% 
      filter(gene %in% gene_subset) %>% 
      filter(P.Value < pvalfilter) %>% 
      filter(logFC > logfcfilter) %>%
      select(gene, logFC, celltype) 
  } )
  mtx = do.call(rbind, mtx)
  mtx = mtx %>% spread(key = celltype, value = logFC) %>% column_to_rownames("gene")
  mtx[is.na(mtx)] = 0
  return(mtx)
}  
## Usage 
#mtx = GetGeneMatrix(result.list = d1time_res, stat_for_matrix = "logFC",gene_subset = core1_genes,pvalfilter = 0.05, logfcfilter = 0)
# pheatmap::pheatmap(mtx)


# make a table of top DE genes 
MakeGeneTable = function(result.list, celltypes.vector, pvalthresh, logfcgreater) { 
  merged = data.frame()
  for (i in 1:length(result.list)) {
    sub = 
      result.list[[i]] %>% 
      mutate(celltype = celltypes.vector[i]) %>% 
      filter(P.Value  < pvalthresh) %>% 
      filter(logFC > logfcgreater)
    merged = rbind(merged, sub)
  }
  return(merged)
}



# ##################################  GENE SET ENRICHMENT ##################################  
# if not using lme4 (e.g. for baseline delta) Get the rank results for the ebayes model fit. 
GetRankResults = function(limma.fit.object.list, coefficient.number, contrast.name){
  print(" returning ranks based on emperical bayes moderated t statistic")
  require(limma)
  ranks = list()
  for (i in 1:length(limma.fit.object.list)) {
    test = as.data.frame(topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf))
    ranks[[i]] = 
      test %>% rownames_to_column("gene") %>%  
      arrange(desc(t)) %>% 
      select(gene, t) %>% 
      column_to_rownames("gene") %>% 
      t() %>% 
      unlist(use.names = T)
    ranks[[i]] = ranks[[i]][1, ]
  }
  names(ranks) = names(limma.fit.object.list)
  return(ranks)
}


GetRankResultsRaw = function(contrast.result.raw.list){
  print("returning genes ranked by t statistic for each cell type based on mixed model results")
  ranks = list()
  for (i in 1:length(contrast.result.raw.list)) {
    ranks[[i]] = contrast.result.raw.list[[i]] %>% 
      arrange(desc(t)) %>% 
      select(gene, t) %>% 
      column_to_rownames("gene") %>% 
      t() %>% 
      unlist(use.names = T)
    ranks[[i]] = ranks[[i]][1, ]
  }
  names(ranks) = names(contrast.result.raw.list)
  return(ranks)
}

# usage: 
# d1hvl_rank = get.rank.results(limma.fit.object.list  = ebf, coefficient.number = 1, contrast.name = "day_1_highvslow")
# d1time_rank = get.rank.results(limma.fit.object.list = ebf, coefficient.number = 3, contrast.name = "day_1_timepoint")

# wrapper around fgsea function 
RunFgseaOnRankList = function(rank.list.celltype, pathways = btm, celltypes.vector = celltypes,
                              maxSize = 500, minSize = 9, nperm = 250000, positive.enrich.only = FALSE) {
  require(fgsea)
  gsea = list()
  for (i in 1:length(rank.list.celltype)) {
    gsea[[i]] = 
      suppressMessages(fgsea(pathways = pathways,
                             stats = rank.list.celltype[[i]],
                             maxSize = maxSize, 
                             minSize = minSize, 
                             nperm = nperm)) %>% 
      mutate(celltype = celltypes.vector[i]) %>%
      arrange(pval)
    if (positive.enrich.only == TRUE) {
      gsea[[i]] = gsea[[i]] %>% filter(ES > 0)
    }
  }
  names(gsea) = celltypes.vector
  return(gsea)
}


## #For plotting fGSEA results in a heatmap, this is copied from Yuri, can also use bubble plot next section to not plot non significant enrichments. 
# MergeGseaResultListFDR = function(gsea.result.list, contrast.name, FDR.threshold, positive.enrich.only = TRUE) { 
#   res.all = data.frame()
#   for (i in 1:length(gsea.result.list)) {
#     if (positive.enrich.only == TRUE) {
#       df = gsea.result.list[[i]] %>% 
#         filter(NES > 0) %>% 
#         select(pathway, padj, celltype) %>%
#         mutate(contrast = rep(contrast.name))
#       res.all = rbind(res.all, df)
#     } else { 
#       df = gsea.result.list[[i]] %>%
#         select(pathway, padj, celltype) %>%
#         mutate(contrast = rep(contrast.name))
#       res.all = rbind(res.all, df)
#     }
#   }
#   threshold = FDR.threshold
#   df.pw =
#     res.all %>% 
#     group_by(pathway) %>%
#     summarise(min.p = min(padj, na.rm=T)) %>%
#     ungroup() %>% 
#     dplyr::filter(min.p < threshold)
#   df.res = 
#     res.all %>%
#     filter(pathway %in% df.pw$pathway)
#   pl = df.res %>% spread(key = celltype, value = padj)
#   pl = pl %>% select(-contrast)
#   plot = pl %>% column_to_rownames("pathway") %>% as.matrix()
#   plot = -log10(plot)
#   plot[is.na(plot)] = 0
#   return(plot)
# }


# MergeGseaResultListNES = function(gsea.result.list, contrast.name, FDR.threshold) { 
#   
#   df = lapply(gsea.result.list, function(x){ x = x %>% select(pathway, padj, NES, celltype) %>% mutate(contrast = rep(contrast.name)) })
#   df = do.call(rbind, df)
#   df = df %>% filter(padj < FDR.threshold)
#   df  = df %>% select(pathway, celltype, NES)
#   
#   pl = df %>% spread(key = celltype, value = NES) %>% column_to_rownames("pathway") %>% as.matrix()
#   # plot = -log10(plot)
#   pl[is.na(pl)] = 0
#   return(pl)
# }



# For GSEA bubble plot instead of heatmap (plots NES and P value and only + enrichment shown)
# merge GSEA list into dataframe 
RbindGseaResultList = function(gsea_result_list, NES_filter = 0, padj_filter = 0.05){
  
  score = lapply(gsea_result_list, function(x){ 
    x = x %>% 
      select(pathway, padj, NES, celltype) %>% 
      filter( NES > NES_filter ) 
  })
  score = do.call(rbind, score) %>% filter(padj < padj_filter) %>% mutate(n_logp = -log10(padj))
  return(score)
}




GSEABubblePlot2 = function(rbind_gsea_result_dataframe, save_path, save_name, width = 8.5, height = 7.2) { 
  p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = n_logp, size = NES)) + 
    geom_point(shape = 21) + 
    scale_fill_viridis_c() + 
    theme_bw() +
    scale_x_discrete(position = "top") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "right") + 
    labs(fill = '-log10(ajdusted P value)', size = 'Normalized Enrichment Score') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
  print(p)
}

GSEABarPlot = function(rbind_gsea_result_dataframe, celltype_name, save_path, title, save_name, fill_color, width = 8.5, height = 7.2) { 
  dplot = rbind_gsea_result_dataframe %>%   
    filter(celltype == celltype_name)
  
  p = ggplot(dplot, aes(y = NES, x = pathway)) + 
    geom_bar(stat="identity", fill=fill_color, alpha=0.8, width= 0.75) +
    coord_flip() +
    theme_bw() +
    scale_fill_viridis_c(option = "B") +
    scale_x_discrete(position = "top") + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.position = "bottom") +  
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) + 
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) +  
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + 
    ggtitle(title)
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
    print(p)
}


GSEABubblePlot = function(rbind_gsea_result_dataframe, save_path, include_negative = FALSE, save_name, width = 8.5, height = 7.2) { 
  
  if (include_negative==TRUE) {
    p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = NES, size = n_logp)) + 
      geom_point(shape = 21) + 
      scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0) +
      theme_bw() +
      scale_x_discrete(position = "top") + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
      theme(axis.title.y = element_blank()) +
      # theme(legend.position = "bottom") + 
      labs(fill = 'Normalized \n Enrichment \n Score', size = '-log10(FDR)') +
      theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
      theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
      theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
      guides(shape = guide_legend(override.aes = list(size = 5))) + 
      guides(color = guide_legend(override.aes = list(size = 5))) + 
      
      theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
    print(p)
  } else{
  p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = n_logp, size = NES)) + 
    geom_point(shape = 21) + 
    scale_fill_viridis_c() + 
    theme_bw() +
    scale_x_discrete(position = "top") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "bottom") + 
    labs(fill = '-log10(ajdusted P value)', size = 'Normalized Enrichment Score') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  }
}



GSEABubblePlot_t = function(rbind_gsea_result_dataframe, save_path, save_name, width = 8.5, height = 7.2) { 
  p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = n_logp, size = NES)) + 
    geom_point(shape = 21) + 
    scale_fill_viridis_c() + 
    theme_bw() +
    scale_x_discrete(position = "top") + 
    # theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "bottom") + 
    labs(fill = '-log10(ajdusted P value)', size = 'Normalized Enrichment Score') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    coord_flip()
  ggsave(p, filename = paste0(save_path,save_name,".png"), width = width, height = height)
  print(p)
  
}

# usage 
# score = RbindGseaResultList(gsea_result_list = gsea1,NES_filter = 0,padj_filter = 0.05)
# GSEABubblePlot2(rbind_gsea_result_dataframe = score,save_path = figpath, save_name = "/btm_0.05",width = 7, height = 5) 


######## COMBINE GSEA AND DE MODEL RESULTS into a single dataframe by gene/ module by p value cutoffs. 
CombineResults = function(gsealist, contrastlist, gseafdr, genefdr){ 
  
  #### 
  '%ni%' = Negate('%in%')
  
  # filter gsea esults and cat pathway and celltype 
  gs = lapply(gsealist, function(x){ x = x %>% filter(padj < gseafdr) %>% mutate(name = paste(pathway, celltype, sep = "__"))})
  
  # get a vector of celltypes with DE results 
  generes = do.call(rbind, contrastlist)
  generes = generes %>% filter(adj.P.Val < genefdr)
  celltype_result = generes$celltype %>% unique()
  
  # subset GSEA results by cell types with genes DE < fdr threshold  
  gs = gs[celltype_result]
  
  # remove any celltypes without gsea results passing fdr 0.01 filter 
  enriched_dim = sapply(gs, dim)
  remove = which(enriched_dim[1, ] == 0) %>% names 
  
  # subset the DE gene list by the celltypes with BTM results 
  generes = generes %>% filter(celltype %ni% remove)
  
  # subset the gsea list by the celltypes with DE genes.
  celltype_result = celltype_result[celltype_result %ni% remove]
  gs = gs[celltype_result]
  
  
  #combine data gs and generes have the same celltype info 
  stopifnot(names(gs) == unique(generes$celltype))
  resdf = list()
  for (i in 1:length(gs)) {
    # get data from 1 cell type 
    #i = 1 
    enrichtest = gs[[i]]  
    gene_celltype = generes %>% filter(celltype == names(gs)[i]) 
    
    # intersect gsea leading edge genes with limma model results 
    lst = apply(enrichtest, 1, function(x) { 
      y = x$leadingEdge %>% as.vector 
      z = gene_celltype %>% 
        filter(gene %in% y) %>% 
        mutate(pathway = rep(x$pathway)) %>% 
        select(gene, celltype, pathway, logFC, padj = adj.P.Val)
      return(z)
    })
    merged = do.call(rbind, lst)
    resdf[[i]] = merged
  }
  result.dataframe = do.call(rbind, resdf) %>% group_by(celltype)
  return(result.dataframe)
}


### Usage: 
#d1time_res = readRDS("mid_res/1_H1N1_pseudobulk_DE/data/d1time_results.rds")
#testgsea = lapply(gsea[[2]], function(x){ x = x %>% filter(NES >0)})
#testmod = lapply(d1time_res, function(x){x %>% filter(logFC > 0)})
#test = CombineResults(gsealist = testgsea, contrastlist = testmod, gseafdr = 0.05,genefdr = 0.1)



# Get the average of leading edge genes from an average data frame list and gsea result list. 
#  most general version for projects. 
GetLeAvg = function(av.exprs.list, gsea.list, padj.filter, NES.filter){ 
  gseasub = lapply(gsea.list, function(x){x = x %>% 
    filter(NES > NES.filter & padj < padj.filter) %>%
    select(pathway, leadingEdge) 
  })
  # get lists with results passing filter and subset result lists by that index  
  g =  which(map_int(gseasub, nrow) > 0) %>% names
  gseasub = gseasub[g] ; avsub = av.exprs.list[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)
  
  mods = lapply(gseasub, function(u){
    testcase = u; paths = u$pathway; dataret = data.frame()
    for (i in 1:length(paths)) {
      genes_use = testcase %>% filter(pathway == paths[i]) %$% leadingEdge %>% unlist %>% as.character()
      dat =  tibble(gene = genes_use, module = rep(paths[i]))
      dataret = rbind(dataret, dat)
    }
    return(dataret)
  })
  tidy = list()
  for (i in 1:length(mods)) {
    average_data = avsub[[i]]
    
    index1 = colnames(average_data)[1]; index2 = colnames(average_data)[ncol(average_data)]
    tidy[[i]] = average_data %>% rownames_to_column("gene") %>% filter(gene %in% mods[[i]]$gene) 
    tidy[[i]] = full_join(tidy[[i]], mods[[i]], by = "gene") %>%
      select(module, everything()) %>% 
      gather(sample, av_exp, index1:index2)
    tidy[[i]]$celltype = names(avsub[i])
  }
  names(tidy) = names(avsub)
  return(tidy)
}


# usage. 
#gsea1 = readRDS("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/data/joint_d1btm.rds")
#av = readRDS("mid_res/1a_Aggregate/data/average_expression_list.rds")
#av = lapply(av, function(x) x %>% select(-contains("d7")))
#tester = GetLeAvg(av.exprs.list = av,gsea.list = gsea1, padj.filter = 0.05, NES.filter = 0)

##### Example customization for plotting: 
# high.responders = c("205","207","209","212","215","234","237","245","250","256")
# test2 = lapply(tester, function(x){
#   x = x %>% 
#     rename("sx_ct" = sample) %>% 
#     mutate(sample = str_sub(sx_ct, -6, -1)) %>% 
#     mutate(cohort = if_else(str_sub(sample, -11, -8) == "H5N1", true =  "H5N1", false = "H1N1")) %>%
#     mutate(timepoint = str_sub(sample, -2,-1)) %>%
#     mutate(group = 
#              if_else(cohort == "H5N1", true = "H5 Adjuvant", 
#                      if_else(cohort == "H1N1" & str_sub(sample, -6,-4) %in% high.responders, true = "H1 high responder", false = "H1 low responder")))  
#   #x$group = factor(x$group, levels = c("H1 low responder" , "H1 high responder" ,"H5 Adjuvant"))
#   
# })

##### 




# get average of top DE genes as tidy data frame 
GetTopAvg = function(av.exprs.list, result.list, P.Value.filter, logFC.filter, top_n_genes = NULL){ 
  
  resultsub = lapply(result.list, function(x){
    x = x %>% 
      filter(logFC > logFC.filter & P.Value < P.Value.filter) %>% 
      arrange(logFC) %$%
      gene 
  })
  
  
  # # get top n genes if specified. 
  if (!is_null(top_n_genes)){
    resultsub = lapply(resultsub, function(x){ x[1:top_n_genes] })
  }else{
    resultsub = resultsub
  }
  
  # get lists with results passing filter and subset result lists by that index  
  g =  which(map_int(resultsub, length) > 0) %>% names
  resultsub = resultsub[g] ; avsub = av.exprs.list[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)
  
  # init list 
  tidy = list()
  for (i in 1:length(avsub)) {
    gene_subset = resultsub[[i]]
    average_data = avsub[[i]][ gene_subset,  ]
    
    index1 = colnames(average_data)[1]; index2 = colnames(average_data)[ncol(average_data)]
    tidy[[i]] = average_data %>% 
      rownames_to_column("gene") %>% 
      # filter(gene %in% resultsub[[i]]) %>% 
      gather(sample, av_exp, index1:index2)
    
    tidy[[i]]$celltype = names(avsub[i])
  }
  names(tidy) = names(avsub)
  return(tidy)
}

## usage: 
# topav = GetTopAvg(av.exprs.list = av,result.list = d1time, P.Value.filter = 0.01,logFC.filter = 1.2)
# 
# # customization of the general results above for this project 
# high.responders = c("205","207","209","212","215","234","237","245","250","256")
# topav = lapply(topav, function(x){
#   x = x %>% 
#     mutate(cohort = if_else(str_sub(sample, -11, -8) == "H5N1", true =  "H5N1", false = "H1N1")) %>%
#     mutate(timepoint = str_sub(sample, -2,-1)) %>%
#     mutate(group = 
#              if_else(cohort == "H5N1", true = "H5 Adjuvant", 
#                      if_else(cohort == "H1N1" & str_sub(sample, -6,-4) %in% high.responders, true = "H1 high responder", false = "H1 low responder")))  
#   x$group = factor(x$group, levels = c("H1 low responder" , "H1 high responder" ,"H5 Adjuvant"))
#   return(x)
# })



#### To do the above on 1 module 1 celltype at at time (for more control over plote etc. )
## LeadingEdge plots 
GetLeadingEdgeGenes = function(gsea.result.list, celltype.index, module.name) {
  
  celltype = gsea.result.list[[celltype.index]]
  print("vector of leadingEdge genes for : ")
  print(unique(celltype$celltype))
  print(module.name)
  genes_use = 
    gsea.result.list[[celltype.index]] %>% 
    filter(pathway %in% module.name) %$% 
    leadingEdge %>% 
    unlist %>% 
    as.character()
  return(genes_use)
}



##### HYPERGEOMETRIC TEST 
RunHypergeometricTest = function(result_list, TERM2GENE_dataframe, pval_threshold = 0.05,
                                 logFC_threshold = 0.5, usefdr_threshold = FALSE){
  
  print("input format: list of dataframes (e.g. per cell type) ")
  print("dataframe column names: gene, logFC, adj.P.Val, P.Value ")
  print( "use PlotHypergeometric() to visualize results for all subsets " )
  
  entrez_subset = lapply(result_list, function(x){
    if (isTRUE(usefdr_threshold)) {
      y = x %>% filter(adj.P.Val < pval_threshold  & logFC > logFC_threshold)
    } else { 
      y = x %>% filter(P.Value < pval_threshold  & logFC > logFC_threshold)   
    }
    map = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                keys = y$gene, columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
    dup = duplicated(map$SYMBOL)
    map = map[!dup, ]
    map = map$ENTREZID
    return(map)
  })
  # remove any unmapped features 
  entrez_subset = lapply(entrez_subset, function(x) x = x[!is.na(x)])
  #init store 
  hypergeometric = enrls =  list()
  for (i in 1:length(entrez_subset)){ 
    
    cells = names(entrez_subset[i])
    print( paste0("hypergeometric test in: ", cells, " index ", i) )
    
    # return NA if here were no genes passing filter above 
    if(length(entrez_subset[[i]]) == 0){ 
      hypergeometric[[i]] = NA 
    } else {  
      # run enrichment 
      enrls[[i]] = 
        suppressMessages(tryCatch(
          clusterProfiler::enricher(entrez_subset[[i]], TERM2GENE = TERM2GENE_dataframe),
          error = function(e) return(NA)
        )) 
      if (is.null(enrls[[i]])) {
        # return NA if there were no enriched pathways within the features in subset i 
        hypergeometric[[i]] = NA
      } else { 
        # reformat enrichment results as dataframe 
        hypergeometric[[i]] = enrls[[i]]
        hypergeometric[[i]] = hypergeometric[[i]]@result %>% 
          mutate(celltype = cells) %>% 
          separate(GeneRatio,into = c("gene_num", "gene_denom"), sep = "/") %>% 
          mutate(gene_num = as.numeric(gene_num)) %>% 
          mutate(gene_denom = as.numeric(gene_denom)) %>% 
          mutate(gene_ratio = round(gene_num / gene_denom, 3))
      }
    }
  }
  # combine results and format for PlotHypergeometric() function
  hyp = do.call(rbind, hypergeometric)
  hyp = hyp %>% mutate(n_logp = -log10(p.adjust))
  return(hyp)
}

# plot hypergeometric results 
PlotHypergeometric = function(hyperg_result, p.adjust.filter = 0.05, genenumber_filter = 0, savepath = figpath, savename , title, height = 10, width = 8 ){ 
  # plot attributes 
  hyperg_attr = list(
    geom_point(shape = 21),
    scale_fill_viridis_c(option = "B"),
    theme_bw(),
    scale_x_discrete(position = "top"),
    theme(axis.text.x=element_text(angle = 45, hjust = 0)),
    theme(axis.title.y = element_blank()),
    # theme(legend.position = "bottom"),
    labs(fill = '-log10(ajdusted P value)', size = 'gene ratio'),
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)),
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")),
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")),
    guides(shape = guide_legend(override.aes = list(size = 5))),
    guides(color = guide_legend(override.aes = list(size = 5))),
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + 
    theme(plot.title = element_text(face = "bold", color = "black", size = 8))
  )
  hyp = hyperg_result 
  p = ggplot(hyp %>% 
               filter(gene_num > genenumber_filter) %>% 
               filter(p.adjust < p.adjust.filter),
             aes(y = ID, x = celltype, fill = n_logp, size = gene_ratio)) + 
    hyperg_attr + 
    ggtitle(title)
  ggsave(p, filename = paste0(savepath, savename, ".pdf"), width = width, height = height)
}



# just run a hypergeometric teste (with the btm term 2 df list e.g.) on a vector of genes - simple 


##### plot average gene distributions for each sample in each cohort. 
# most general version 
GetTidySummary = function(av.exprs.list, celltype.index, genes.use, subset_cohort = FALSE, cohort.plot = NULL){ 
  gene.index.1 = genes.use[1]
  gene.index.2 = genes.use[length(genes.use)]
  tidy_data = 
    av.exprs.list[[celltype.index]][genes.use, ] %>% 
    t %>% 
    as.data.frame() %>% 
    rownames_to_column("sample") %>% 
    gather(key = gene, value = count, gene.index.1:gene.index.2)
  return(tidy_data)
}



# plot the tidy gene x sample summary by cohort. 
PlotGeneDistCohort = function(merged_av_data, save_path, save_name, title = NULL, nrow = 5, height = 8, width = 10, plot_subset = FALSE, genes_plot = NULL){
  print("merged av data needs to have a group column and a timepoint column if color error add more colors to cu = in first line of function")
  #cu = pals::stepped()
  # cu = cu[c(12,9,16,13,4,1)]
  #cu = cu[c(12,9,16,13)]
  # cu = c(cu,  "#FFBB78FF", "#FF7F0EFF")
  cu = c("dodgerblue", "midnightblue", "red", "firebrick",  "#FFBB78FF", "#FF7F0EFF")
  
  if(plot_subset == TRUE) {
    p = ggplot(merged_av_data %>% filter(gene %in% genes_plot), aes(x = group, y = count, fill = interaction(timepoint, group ))) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free_y", nrow = nrow) +
      theme_bw(base_size = 10.5)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      # scale_fill_brewer(palette = 'Paired') +
      scale_fill_manual(values = cu) +
      # theme(strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid")) +
      theme(strip.background = element_blank()) +
      theme(strip.text = element_text(face = "bold",family = "Helvetica")) +
      theme(axis.text.x =  element_blank()) +
      theme(axis.text.y =  element_text(size = 6)) +
      ggtitle(title)+ 
      ylab("Average Expression")
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), height = height,  width = width)
    
  } else {
    p = ggplot(merged_av_data, aes(x = group, y = count, fill = interaction(timepoint, group ))) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free_y", nrow = nrow) +
      theme_bw(base_size = 10.5)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      # scale_fill_viridis_d() +
      # scale_fill_brewer(palette = 'Paired') +
      scale_fill_manual(values = cu) +
      theme(strip.background = element_blank()) +
      # theme(strip.background = element_rect(color="black", fill="white", size=0.5, linetype="solid")) +
      theme(strip.text = element_text(face = "bold",family = "Helvetica")) +
      theme(axis.text.x =  element_blank()) +
      theme(axis.text.y =  element_text(size = 6)) +
      ggtitle(title)+ 
      ylab("Average Expression")
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), height = height,  width = width)
  }
}


## Just Get the leading edge genes for the full dataset. 
GetLeadingEdgeFull = function(gsea.list, padj.filter, NES.filter){ 
  gseasub = lapply(gsea.list, function(x){x = x %>% 
    filter(NES > NES.filter & padj < padj.filter) %>%
    select(pathway, leadingEdge) 
  })
  # get lists with results passing filter and subset result lists by that index  
  g =  which(map_int(gseasub, nrow) > 0) %>% names
  gseasub = gseasub[g] 
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)
  mods = lapply(gseasub, function(u){
    testcase = u; paths = u$pathway; dataret = data.frame()
    for (i in 1:length(paths)) {
      genes_use = testcase %>% filter(pathway == paths[i]) %$% leadingEdge %>% unlist %>% as.character()
      dat =  tibble(gene = genes_use, module = rep(paths[i]))
      dataret = rbind(dataret, dat)
    }
    return(dataret)
  })
  celltypes = names(gseasub)
  for (i in 1:length(mods)) {
    mods[[i]]$celltype = celltypes[i]
    }
  return(mods)
}
# usage: 
# add wrapper to GetTidySummary for plots to specify timepoint cohort, response group.
# GetTidyCohort = function(av.exprs.list , celltype.index , genes.use ){
#   plot_df = 
#     GetTidySummary(av.exprs.list = av.exprs.list, celltype.index = celltype.index, genes.use = genes.use) %>%
#     mutate(cohort = if_else(str_sub(sample, -11, -8) == "H5N1", true =  "H5N1", false = "H1N1")) %>%
#     mutate(timepoint = str_sub(sample, -2,-1)) %>%
#     mutate(group = 
#              if_else(cohort == "H5N1", true = "H5 Adjuvant", 
#                      if_else(cohort == "H1N1" & str_sub(sample, -6,-4) %in% high.responders, true = "H1 high responder", false = "H1 low responder")))
#   plot_df$group= factor(plot_df$group, levels = c("H1 low responder" , "H1 high responder" ,"H5 Adjuvant"))
#   return(plot_df)
# }
# 
# 
# ## CD14 Monocytes (index 7) TLR BTM
# mono_tlr = GetLeadingEdgeGenes(gsea.result.list = gsea1, celltype.index = 7, module.name = gsea1[[7]][3,]$pathway)
# mtlr = GetTidyCohort(av.exprs.list = av,celltype.index = 7, genes.use = mono_tlr)
# PlotGeneDistCohort(merged_av_data = mtlr, save_name = "monocyte_TLR", save_path = figpath,  height = 4.7 , width = 16, nrow = 2)
# 


### match order of 2 lists of data frames based on "celltype" column to match fgsea result lists 
MatchfgseaResultIndex = function(list_to_reference, list_to_reorder){
  
  print("cell type names must match")
  
  result.list1 = list_to_reference
  result.list2 = list_to_reorder
  
  ct = sapply(result.list1, function(x){ unique(x$celltype)}) %>% unname
  print("list 1 order of results by cell type")
  print(ct)
  ct2 = sapply(result.list2, function(x){ unique(x$celltype) }) %>% unname
  print("list 2 order of results by cell type")
  print(ct2)
  
  if (!(length(ct) == length(ct2))) {
    print("the celltype result lists are unequal length and will be merged by the intersection")
  }
  
  ct2 = ct2[ct2 %in% ct]
  ct = ct[ct %in% ct2]
  
  reorder_2_match1 = ct2[order(match(ct2, ct))]
  
  result.list2 = result.list2[reorder_2_match1]
  return(result.list2)
}



# Make volcano plot for each celltype 
VolcanoPlotTop = function(contrast.result.list, contrast.name, save.path, size = 3, fig.height, fig.width) {
  require(ggrepel)
  for (i in 1:length(contrast.result.list)) {
    #format result list for volcano plot 
    contrast.result.list[[i]] = contrast.result.list[[i]] %>% 
      arrange(adj.P.Val)  %>% 
      mutate(nlp = -log10(adj.P.Val)) %>% 
      mutate(color = 
               if_else(nlp > 1 & LogFC > 0.2, true = "1", 
                       if_else(nlp > 1 & LogFC < -0.2,  true = "2", 
                               if_else(nlp > 1 & LogFC < 0.2 & LogFC > - 0.2,  true = "3", false = "0"))))
    
    # volcano plot 
    p = ggplot(data = contrast.result.list[[i]], aes(x = LogFC, y = nlp, label = gene, color = color)) + 
      geom_point(show.legend = FALSE) + 
      theme_bw() +
      geom_text_repel(data= contrast.result.list[[i]] %>% 
                        filter(nlp > 1.3) %>% 
                        filter(abs(LogFC) > 0.2),  
                      aes(label=gene, color = color), size = size ,show.legend = FALSE) +
      geom_vline(xintercept = 0.2, linetype = "dashed") + 
      geom_vline(xintercept =  - 0.2, linetype = "dashed") + 
      geom_hline(yintercept = 1.3, linetype = "dashed") + 
      labs(x = "Log FC", y = "-log10 adjusted p value") +
      scale_color_manual(values = c("black", "red", "blue", "black")) + 
      theme(text =  element_text(colour = "black", family = "Helvetica",size = 12)) + 
      ggtitle(names(contrast.result.list[i]))
    
    # save plot 
    ggsave(plot = p, 
           filename = paste0(contrast.name, names(contrast.result.list)[i],".pdf"),
           path = save.path, 
           width = fig.width, height = fig.height)
  }
}


