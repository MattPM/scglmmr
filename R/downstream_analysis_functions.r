#### De functions not dependent on Seurat versions. Downstream enrichment etc.
# select = dplyr::select
'%ni%' = Negate('%in%')

#' GetContrastResults - return results from a contrast fit on list of celltypes from RunVoomLimma using topTable
#'
#' @param limma.fit.object.list the results returned by RunVoomLimma, to get coefficiennt from dreamMixedModel use GetContrastResultsRaw
#' @param coefficient.number corresponds to the contrast, the nmber is in order of the contrast matrix
#' @param contrast.name this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.
#' @return
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate if_else
#' @export
#'
#'
#' @examples
GetContrastResults = function(limma.fit.object.list, coefficient.number, contrast.name){
  print("this function returns results from RunVoomLimma, to get coefficient from dreamMixedModel, use GetContrastResultsRaw
        GetContrastResults uses emperican Bayes shrinkage see https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  # init store
  test = list()

  for (i in 1:length(limma.fit.object.list)) {
    test[[i]] =
      limma::topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf, sort.by = "p") %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(contrast = rep(contrast.name)) %>%
      dplyr::mutate(pvalue0.01 = dplyr::if_else(P.Value < 0.01, true = "1", false = "0")) %>%
      dplyr::mutate(celltype = celltypes.vector[i])
  }
  names(test) = names(limma.fit.object.list)
  return(test)
}

#' GetContrastResultsRaw - calculate p values and return contrast results from modelfit with dreamMixedModel
#'
#' @param limma.fit.object.list the results returned by dreamMixedModel, to get coefficiennt from RunVoomLimma use GetContrastResults
#' @param coefficient.number corresponds to the contrast, the nmber is in order of the contrast matrix
#' @param contrast.name this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.
#' @return
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate full_join select
#' @importFrom stats p.adjust
#' @importFrom plyr mapvalues
#' @export
#'
#' @examples
GetContrastResultsRaw =  function(limma.fit.object.list, coefficient.number, contrast.name){
  print("this function returns results from dreamMixedModel, to get coefficient from RunVoomLimma, use GetContrastResults
        GetContrastResults uses emperican Bayes shrinkage see https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  ## pt 1 return ONLY gene, logFC, AveExpr, contrast, celltype from eBayes call to format data and add raw p values from lme4, correct with BH.
  pvals = lapply(limma.fit.object.list, function(x){ x$pValue[ , coefficient.number] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::rename(P.Value = '.') })

  # parameters to run ordinary t statistic
  coef = lapply(limma.fit.object.list, function(x){ x$coefficients[ ,coefficient.number] })
  stdev = lapply(limma.fit.object.list, function(x){ x$stdev.unscaled[ ,coefficient.number] })
  sigma_ = lapply(limma.fit.object.list, function(x){ x$sigma })
  # Note eBayes t statistics are NOT used, only for formatting output, next section adds raw p values and logFC returned by lme4
  ebf = lapply(limma.fit.object.list, eBayes) ##
  contrast_result = GetContrastResults(limma.fit.object.list = ebf,
                                       coefficient.number = coefficient.number,
                                       sort.by = "p",
                                       contrast.name = contrast.name)
  contrast_result = lapply(contrast_result, function(x){ x %>% dplyr::select(gene, logFC, AveExpr, contrast, celltype) })
  result = t_statistic = list()
  for (i in 1:length(limma.fit.object.list)) {
    # Ordinary t-statistic (see e.g. limma eBayes documentation: ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
    t_statistic[[i]] = coef[[i]] / stdev[[i]] / sigma_[[i]] %>% as.data.frame()
    t_statistic[[i]] =  t_statistic[[i]] %>% tibble::rownames_to_column("gene")
    colnames(t_statistic[[i]]) = c("gene", "t")

    result[[i]] =
      contrast_result[[i]] %>%
      dplyr::mutate(P.Value = plyr::mapvalues(x = gene, from = pvals[[i]]$gene, to = round(pvals[[i]]$P.Value, digits = 9))) %>%
      dplyr::mutate(adj.P.Val = stats::p.adjust(p = P.Value, method = "BH"))
    # add t statistic
    result[[i]] = dplyr::full_join(result[[i]], t_statistic[[i]], by = "gene")

  }
  names(result) = names(limma.fit.object.list)
  return(result)
}


#' GetGeneMatrix get a gene matric for plotting genes by celltypes statistic from pseudobulk model results e.g. using heatmap pheatmap or complexheatmap
#'
#' @param result.list the object returned by GetContrastResults or GetContrastResultsRaw
#' @param gene_subset a preselected subset of genes as a vector, for example returned by GetLeadingEdgeFull. Defaults to the union of all fitted genes across cell types.
#' @param stat_for_matrix defaults to logFC, the effect size. can be any of the columns returned by limma::topTable
#' @param pvalfilter filter genes to retain in the matrix by the raw p values
#' @param logfcfilter filter genes to retain in the matrix by the logFC
#'
#' @return
#' @importFrom dplyr filter select bind_rows
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
GetGeneMatrix = function(result.list, gene_subset, stat_for_matrix = "logFC", pvalfilter, logfcfilter){
  gene_subset = do.call(rbind, test)$gene %>% unique
  mtx = lapply(result.list, function(x){
    x = x %>%
      dplyr::filter(gene %in% gene_subset) %>%
      dplyr::filter(P.Value < pvalfilter) %>%
      dplyr::filter(logFC > logfcfilter) %>%
      dplyr::select(gene, logFC, celltype) }) %>%
      dplyr::bind_rows() %>%
      tidyr::spread(key = celltype, value = logFC) %>%
      tibble::column_to_rownames("gene")
  mtx[is.na(mtx)] = 0
  print("NA values for cell types that did not have a model fit for a given gene are converted to 0")
  return(mtx)
}

# ##################################  GENE SET ENRICHMENT ##################################


#' GetRankResults get list of gene ranks by t stat for fGSEA
#'
#' @param limma.fit.object.list the results returned by RunLimmaVoom. Use GetRankResultsRaw for results returned by dreamMixedModel
#' @param coefficient.number
#' @param contrast.name
#'
#' @return
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr arrange select
#' @importFrom limma topTable
#' @export
#'
#' @examples
GetRankResults = function(limma.fit.object.list, coefficient.number, contrast.name){
  print(" returning ranks based on emperical bayes moderated t statistic")
  ranks = list()
  for (i in 1:length(limma.fit.object.list)) {
    test = as.data.frame(limma::topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf))
    ranks[[i]] =
      test %>% tibble::rownames_to_column("gene") %>%
      dplyr::arrange(dplyr::desc(t)) %>%
      dplyr::select(gene, t) %>%
      tibble::column_to_rownames("gene") %>%
      t() %>%
      unlist(use.names = T)
    ranks[[i]] = ranks[[i]][1, ]
  }
  names(ranks) = names(limma.fit.object.list)
  return(ranks)
}


#' GetRankResults get list of gene ranks by t stat for fGSEA. for results returned by dreamMixedModel
#'
#' @param contrast.result.raw.list results returned by dreamMixedModel
#'
#' @return
#' @importFrom dplyr arrange select desc
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
GetRankResultsRaw = function(contrast.result.raw.list){
  print("returning genes ranked by t statistic for each cell type based on mixed model results")
  ranks = list()
  for (i in 1:length(contrast.result.raw.list)) {
    ranks[[i]] = contrast.result.raw.list[[i]] %>%
      dplyr::arrange(dplyr::desc(t)) %>%
      dplyr::select(gene, t) %>%
      tibble::column_to_rownames("gene") %>%
      t() %>%
      unlist(use.names = T)
    ranks[[i]] = ranks[[i]][1, ]
  }
  names(ranks) = names(contrast.result.raw.list)
  return(ranks)
}


#' RunFgseaOnRankList - wrapper around fast gene set enrichment analysis with the fgsea R package https://bioconductor.org/packages/release/bioc/html/fgsea.html
#'
#' @param rank.list.celltype results returned by GetRankResultsRaw or GetRankResults
#' @param pathways modules / gene sets as a named list each a single vector of unique gene IDS
#' @param maxSize see fgsea package
#' @param minSize see fgsea package
#' @param nperm recommended to keep set at 25000 based on optomization and p value stabilization
#' @param positive.enrich.only include negative enrichments in results? TRUE/FALSE
#'
#' @return
#' @importFrom fgsea fgsea
#' @importFrom dplyr arrange filter
#' @export
#'
#' @examples
#' # usage:
#' t1hvl_rank = GetRankResultsRaw(limma.fit.object.list  = ebf, coefficient.number = 1, contrast.name = "time_1_highvslow")
#' gsea = RunFgseaOnRankList(rank.list.celltype = t1hvl_rank, )
RunFgseaOnRankList = function(rank.list.celltype, pathways, maxSize = 500, minSize = 9,
                              nperm = 250000, positive.enrich.only = FALSE) {
  #require(fgsea)
  gsea = list()
  for (i in 1:length(rank.list.celltype)) {
    gsea[[i]] =
      suppressMessages(fgsea::fgsea(pathways = pathways,
                                    stats = rank.list.celltype[[i]],
                                    maxSize = maxSize,
                                    minSize = minSize,
                                    nperm = nperm)) %>%
      dplyr::mutate(celltype = celltypes.vector[i]) %>%
      dplyr::arrange(pval)
    if (positive.enrich.only == TRUE) {
      gsea[[i]] = gsea[[i]] %>% dplyr::filter(ES > 0)
    }
  }
  names(gsea) = names(rank.list.celltype)
  return(gsea)
}


#' RbindGseaResultList - prepare gsea result list for visualization funcitons GseaBubblePlot or GseaBarPlot
#'
#' @param gsea_result_list
#' @param NES_filter
#' @param padj_filter
#'
#' @return
#' @importFrom dplyr select filter mutate
#' @export
#'
#' @examples
RbindGseaResultList = function(gsea_result_list, NES_filter = -Inf, padj_filter = 0.05){
  score = lapply(gsea_result_list, function(x){
    x = x %>%
      dplyr::select(pathway, padj, NES, celltype) %>%
      dplyr::filter( NES > NES_filter )
  })
  score = do.call(rbind, score) %>% dplyr::filter(padj < padj_filter) %>% dplyr::mutate(n_logp = -log10(padj))
  return(score)
}

#' GSEABarPlot - plot gsea results for a single cell type
#'
#' @param rbind_gsea_result_dataframe
#' @param celltype_name
#' @param save_path
#' @param title
#' @param save_name
#' @param fill_color
#' @param width
#' @param height
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
GSEABarPlot = function(rbind_gsea_result_dataframe, celltype_name, save_path, title, save_name, fill_color, width = 8.5, height = 7.2) {
  # filter to single cell type
  dplot = rbind_gsea_result_dataframe[ rbind_gsea_result_dataframe$celltype == celltype_name, ]

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


#' GSEABubblePlot plot gsea results for all cell types
#'
#' @param rbind_gsea_result_dataframe dataframe returned by RbindGseaResultList
#' @param save_path
#' @param include_negative TRUE/FALSE whether to include negative enrichment in the plot.
#' @param save_name
#' @param width
#' @param height
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
GSEABubblePlot = function(rbind_gsea_result_dataframe, save_path, include_negative = TRUE, save_name, width = 8.5, height = 7.2) {

  # apply same aes for both includeneg with exception
  plot_param = list (
    geom_point(shape = 21),
    theme_bw(),
    scale_x_discrete(position = "top"),
    theme(axis.text.x=element_text(angle = 45, hjust = 0)),
    theme(axis.title.y = element_blank()),
    labs(fill = 'Normalized \n Enrichment \n Score', size = '-log10(FDR)'),
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)),
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")),
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")),
    guides(shape = guide_legend(override.aes = list(size = 5))),
    guides(color = guide_legend(override.aes = list(size = 5))),
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  )

  if (include_negative==TRUE) {
    p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = NES, size = n_logp)) +
      plot_param +
      scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0)
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
    print(p)
  } else{
  p = ggplot(rbind_gsea_result_dataframe, aes(y = pathway, x = celltype, fill = n_logp, size = NES)) +
    plot_param +
    scale_fill_viridis_c()
  ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
  print(p)
  }
}



# usage
# score = RbindGseaResultList(gsea_result_list = gsea1,NES_filter = 0,padj_filter = 0.05)
# GSEABubblePlot2(rbind_gsea_result_dataframe = score,save_path = figpath, save_name = "/btm_0.05",width = 7, height = 5)
## Just Get the leading edge genes for the full dataset.


#' GetLeadingEdgeFull Get a tidy dataframe of ALL Leading Edge Genes from gene set enrichment for all cell types
#'
#' @param gsea.list results returned by RunFgseaOnRankList
#' @param padj.filter filter reslts
#' @param NES.filter filter results
#'
#' @return
#' @importFrom purrr map_int
#' @importFrom dplyr filter select
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' y = GetLeadingEdgeFull(gsea, padj.filter = 0.1, NES.filter = -Inf)
GetLeadingEdgeFull = function(gsea.list, padj.filter, NES.filter){
  gseasub = lapply(gsea.list, function(x){x = x %>%
    dplyr::filter(NES > NES.filter & padj < padj.filter) %>%
    dplyr::select(pathway, leadingEdge)
  })
  # get lists with results passing filter and subset result lists by that index
  g =  which(purrr::map_int(gseasub, nrow) > 0) %>% names
  gseasub = gseasub[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)
  mods = lapply(gseasub, function(u){
    testcase = u; paths = u$pathway; dataret = data.frame()
    for (i in 1:length(paths)) {
      genes_use = testcase %>%
        dplyr::filter(pathway == paths[i]) %$% leadingEdge %>%
        unlist %>%
        as.character()
      dat =  tibble::tibble(gene = genes_use, module = rep(paths[i]))
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


#### To do the above on 1 module 1 celltype at at time (for more control over plote etc. )
## LeadingEdge plots


#' GetLeadingEdgeGenes get the leading edge genes from a single cell type / module combination
#'
#' @param gsea.result.list results from RunFgseaOnRank
#' @param celltype.index celltype number of result list
#' @param module.name name of module
#'
#' @return
#' @importFrom dplyr filter
#' @export
#'
#' @examples
GetLeadingEdgeGenes = function(gsea.result.list, celltype.index, module.name) {

  celltype = gsea.result.list[[celltype.index]]
  print("vector of leadingEdge genes for : ")
  print(unique(celltype$celltype))
  print(module.name)
  genes_use =
    gsea.result.list[[celltype.index]] %>%
    dplyr::filter(pathway %in% module.name) %$%
    leadingEdge %>%
    unlist %>%
    as.character()
  return(genes_use)
}


#' CombineResults - For all cell types merge the GeneSetEnrichment leading edge genes coefficieint and p value from limma
#'
#' @param gsealist
#' @param contrastlist
#' @param gseafdr
#' @param genefdr
#'
#' @return
#' @importFrom dplyr mutate filter select group_by
#' @export
#'
#' @examples
#' test = CombineResults(gsealist = testgsea, contrastlist = testmod, gseafdr = 0.05,genefdr = 0.1)
CombineResults = function(gsealist, contrastlist, gseafdr, genefdr){
  '%ni%' = Negate('%in%')

  # filter gsea esults and cat pathway and celltype
  gs = lapply(gsealist, function(x){
    x = x %>%
      dplyr::filter(padj < gseafdr) %>%
      dplyr::mutate(name = paste(pathway, celltype, sep = "__"))})

  # get a vector of celltypes with DE results
  generes = do.call(rbind, contrastlist)
  generes = generes %>% dplyr::filter(adj.P.Val < genefdr)
  celltype_result = generes$celltype %>% unique()

  # subset GSEA results by cell types with genes DE < fdr threshold
  gs = gs[celltype_result]

  # remove any celltypes without gsea results passing fdr filter
  enriched_dim = sapply(gs, dim)
  remove = which(enriched_dim[1, ] == 0) %>% names

  # subset the DE gene list by the celltypes with BTM results
  generes = generes %>% dplyr::filter(celltype %ni% remove)

  # subset the gsea list by the celltypes with DE genes.
  celltype_result = celltype_result[celltype_result %ni% remove]
  gs = gs[celltype_result]

  #combine data gs and generes have the same celltype info
  stopifnot(names(gs) == unique(generes$celltype))
  resdf = list()
  for (i in 1:length(gs)) {
    # get data from 1 cell type
    enrichtest = gs[[i]]
    gene_celltype = generes %>% dplyr::filter(celltype == names(gs)[i])

    # intersect gsea leading edge genes with limma model results
    lst = apply(enrichtest, 1, function(x) {
      y = x$leadingEdge %>% as.vector
      z = gene_celltype %>%
        dplyr::filter(gene %in% y) %>%
        dplyr::mutate(pathway = rep(x$pathway)) %>%
        dplyr::select(gene, celltype, pathway, logFC, padj = adj.P.Val)
      return(z)
    })
    merged = do.call(rbind, lst)
    resdf[[i]] = merged
  }
  result.dataframe = do.call(rbind, resdf) %>% dplyr::group_by(celltype)
  return(result.dataframe)
}


#' LeadEdgeTidySampleExprs - convert a PseudobulkList into a tidy dataframe for each sample across cell types of the leading edge genes from a gsea list
#'
#' @param av.exprs.list object returned by PseudobulkList summed or average counts
#' @param gsea.list object returned by RunFgseaOnRankList
#' @param padj.filter filter for adjusted p from GSEA
#' @param NES.filter filter for normalized enrichment score from GSEA
#'
#' @return
#' @importFrom purrr map_int
#' @importFrom dplyr filter select
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather everything
#' @export
#'
#' @examples
#' pos_enr = GetLeadingEdgeAverage(av.exprs.list = pb_average, gsea.list = gsea1, padj.filter = 0.1, NES.filter = 0)
LeadEdgeTidySampleExprs = function(av.exprs.list, gsea.list, padj.filter, NES.filter){
  gseasub = lapply(gsea.list, function(x){x = x %>%
    dplyr::filter(NES > NES.filter & padj < padj.filter) %>%
    dplyr::select(pathway, leadingEdge)
  })
  # get lists with results passing filter and subset result lists by that index
  g =  which(purrr::map_int(gseasub, nrow) > 0) %>% names
  gseasub = gseasub[g] ; avsub = av.exprs.list[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)

  mods = lapply(gseasub, function(u){
    testcase = u; paths = u$pathway; dataret = data.frame()
    for (i in 1:length(paths)) {
      genes_use = testcase %>% dplyr::filter(pathway == paths[i]) %$% leadingEdge %>% unlist %>% as.character()
      dat =  tibble::tibble(gene = genes_use, module = rep(paths[i]))
      dataret = rbind(dataret, dat)
    }
    return(dataret)
  })
  # pb data
  tidy = list()
  for (i in 1:length(mods)) {
    average_data = avsub[[i]]

    index1 = colnames(average_data)[1]; index2 = colnames(average_data)[ncol(average_data)]
    tidy[[i]] = average_data %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% mods[[i]]$gene)

    # combine with modules
    tidy[[i]] = dplyr::full_join(tidy[[i]], mods[[i]], by = "gene") %>%
      dplyr::select(module, tidyr::everything()) %>%
      tidyr::gather(sample, av_exp, index1:index2)
    tidy[[i]]$celltype = names(avsub[i])
  }
  names(tidy) = names(avsub)
  return(tidy)
}


#' TopGenesTidySampleExprs convert a PseudobulkList into a tidy dataframe for each sample across cell types of the top differentially expressed genes for a contrast
#'
#' @param av.exprs.list
#' @param result.list
#' @param P.Value.filter
#' @param logFC.filter
#' @param top_n_genes
#'
#' @return
#' @importFrom dplyr filter arrange
#' @importFrom purrr map_int
#' @importFrom tibble rownames_to_column
#' @tidyr gather
#' @export
#'
#' @examples
#' topav = GetTopAvg(av.exprs.list = av,result.list = d1time, P.Value.filter = 0.01,logFC.filter = 1.2)
TopGenesTidySampleExprs = function(av.exprs.list, result.list, P.Value.filter, logFC.filter, top_n_genes = NULL){

  resultsub = lapply(result.list, function(x){
    x = x %>%
      dplyr::filter(logFC > logFC.filter & P.Value < P.Value.filter) %>%
      dplyr::arrange(logFC) %$%
      gene
  })

  # # get top n genes if specified.
  if (!is_null(top_n_genes)){
    resultsub = lapply(resultsub, function(x){ x[1:top_n_genes] })
  }else{
    resultsub = resultsub
  }

  # get lists with results passing filter and subset result lists by that index
  g =  which(purrr::map_int(resultsub, length) > 0) %>% names
  resultsub = resultsub[g] ; avsub = av.exprs.list[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)

  # init list
  tidy = list()
  for (i in 1:length(avsub)) {
    gene_subset = resultsub[[i]]
    average_data = avsub[[i]][ gene_subset,  ]

    index1 = colnames(average_data)[1]; index2 = colnames(average_data)[ncol(average_data)]
    tidy[[i]] = average_data %>%
      tibble::rownames_to_column("gene") %>%
      tidyr::gather(sample, av_exp, index1:index2)

    tidy[[i]]$celltype = names(avsub[i])
  }
  names(tidy) = names(avsub)
  return(tidy)
}

## usage:

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



RunHypergeometricTest = function(result_list, TERM2GENE_dataframe, pval_threshold = 0.05,
                                 logFC_threshold = 0.5, usefdr_threshold = FALSE){
  # result_list = resl[[6]]
  # TERM2GENE_dataframe = term_df_hllmark
  # pval_threshold = 0.05
  # logFC_threshold = 0
  # usefdr_threshold = TRUE
  print("result_list = list of dataframes indexed by celltype \n column names: gene, logFC, adj.P.Val, P.Value; use PlotHypergeometric() on results")
  entrez_subset = list()
  for (i in 1:length(result_list)) {

    if (isTRUE(usefdr_threshold)) {
      result_list[[i]] = result_list[[i]] %>% filter(adj.P.Val < pval_threshold  & logFC > logFC_threshold)
    } else {
      result_list[[i]] = result_list[[i]] %>% filter(P.Value < pval_threshold  & logFC > logFC_threshold)
    }
    # map geneID to entrez ids
    entrez_subset[[i]] =
      tryCatch(
        AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                              keys = result_list[[i]]$gene, columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL"),
        error = function(e) return(NA)
      )
    # if no ids are mapped keep the NA value.
    if (is.na(entrez_subset[i])) {
      entrez_subset[[i]] = NA
    } else {
      dup = duplicated(entrez_subset[[i]]$SYMBOL)
      entrez_subset[[i]] = entrez_subset[[i]][!dup, ]
      entrez_subset[[i]] = entrez_subset[[i]]$ENTREZID
    }
  }
  # remove any unmapped features
  entrez_subset = lapply(entrez_subset, function(x) x = x[!is.na(x)])
  names(entrez_subset) = names(result_list)
  # init store
  hypergeometric = list()
  for (i in 1:length(entrez_subset)){

    cells = names(entrez_subset[i])
    print(paste0("hypergeometric test in: ", cells, " index ", i))
    # at least 3 genes
    if(length(entrez_subset[[i]]) <= 3){
      hypergeometric[[i]] = NA
    } else {

      # run enrichment
      hypergeometric[[i]] = suppressMessages(tryCatch(
        clusterProfiler::enricher(entrez_subset[[i]], TERM2GENE = TERM2GENE_dataframe),
        error = function(e) return(NA)))

      # return NA if there were no enriched pathways within the Entrez IDS in hypergeometric[[i]]
      if (is.null(hypergeometric[[i]]@result)) {
        hypergeometric[[i]] = NA
      } else {
        # reformat enrichment results as dataframe
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
  hyp %<>% mutate(n_logp = -log10(p.adjust))
  return(hyp)
}



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


#############
# retired functions

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




