# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com


#' ExtractResult - convenience function to return statistics for downstream analysis functions such as FgseaList. Returns results from list of dream or lmFit results from those functions natively (use model.fit.list = list(fit)) or from scglmmr::RunVoomLimma and scglmmr::dreamMixedModel
#' @param model.fit.list list of model results indexed by celltypes returned by `scglmmr::dreamMixedModel`,  `scglmmr::RunVoomLimma`, or manually by `lmFit` or `dream`.
#' @param what what to return what = c('statistics', 'gene.t.ranks')[1] defaults to statistics for each cell type, e.g. avg exprs, logFC, t statistic, pval, adj.P.Val etc. If gene.t.ranks, ranks genes based on t statistic and returns a named numeric vector for FgseaList.
#' @param coefficient.number what coefficient to return -- this needs to be one of model.fit.list$coefficients: check the order of the coefficients. Results returned from dream include statistical contrasts and estimated coefficients from the model. If limma::contrasts.fit was used (e.g. if using do_contrast_fit = TRUE in RunVoomLimma), these 'coefficients' are results of the statistical contrast.
#' @param coef.name the name of the estimated coefficient for which results are being returned; if returning results from a statistical contrast e.g. limma::contrasts.fit() this will be the name of the contrast. If returning a model fit with the dream function, can also be contrast specified by from variancePartition::makeContrastsDream() or a fixed effect parameter that was included in the model.
#' @return a list of dataframes with contrast results indexed by cell type or a list of genes ranked b t statistic in format ready for FgseaList.
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select
#' @export
ExtractResult = function(model.fit.list, what = c('statistics', 'gene.t.ranks')[1], coefficient.number, coef.name){
  #init
  celltypes = names(model.fit.list)

  # check user input data
  model_type = model.fit.list[[1]]$method
  message1 = 'raw t statistic reported for unequal degrees of freedom'
  message2 = 'emperical Bayes moderated t statistic reported for model with only fixed effects'
  messageprint = ifelse(model_type == 'lmer', message1, message2)
  print(paste0('returning results of ', model_type, ' model: ', messageprint))
  coefs = colnames(model.fit.list[[1]]$coefficients)
  print('coefficients available from model fit object: ');print(coefs)
  # ensure user is estimating the contrast that
  stopifnot(all.equal(
    as.character(coef.name), as.character(coefs[coefficient.number])
  ))

  # extract results
  test = ret = list()
  for (i in 1:length(model.fit.list)) {
    test[[i]] =
      limma::topTable(fit = model.fit.list[[i]], coef = coefficient.number, number = Inf, sort.by = 't') %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(contrast = rep(coef.name)) %>%
      dplyr::mutate(celltype = celltypes[i]) %>%
      arrange(desc(t))

    if (what == 'gene.t.ranks') {
      ret[[i]] = test[[i]] %>%
        dplyr::select(c('gene', 't')) %>%
        tibble::column_to_rownames("gene") %>%
        t() %>%
        unlist(use.names = T)
      ret[[i]] = ret[[i]][1, ]
    } else{
      ret[[i]] = test[[i]]
    }
  }
  names(ret) = names(model.fit.list)
  return(ret)
}


#' GetContrastResults - return results from a contrast fit on list of celltypes from RunVoomLimma using topTable
#'
#' @param limma.fit.object.list the results returned by RunVoomLimma, to get coefficiennt from dreamMixedModel use GetContrastResultsRaw
#' @param coefficient.number corresponds to the contrast, the nmber is in order of the contrast matrix
#' @param contrast.name this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.
#' @return a list of dataframes with contrast results indexed by cell type
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate if_else
#' @export
#'
#'
#' @examples
#'\dontrun{
#'res = scglmmr::GetContrastResults(limma.fit.object.list = bl, coefficient.number = 1, contrast.name = "test")
#' }
GetContrastResults = function(limma.fit.object.list, coefficient.number, contrast.name){
  print('this function will be deprecated by ExtractResult')

  print("this function returns results from RunVoomLimma, to get coefficient from dreamMixedModel, use GetContrastResultsRaw
        GetContrastResults uses emperican Bayes shrinkage see https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  # init store
  test = list()
  celltypes = names(limma.fit.object.list)

  for (i in 1:length(limma.fit.object.list)) {
    test[[i]] =
      limma::topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf, sort.by = "p") %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(contrast = rep(contrast.name)) %>%
      dplyr::mutate(celltype = celltypes[i])
  }
  names(test) = names(limma.fit.object.list)
  return(test)
}

#' GetContrastResultsRaw - calculate p values and return contrast results from modelfit with dreamMixedModel
#'
#' @param limma.fit.object.list the results returned by dreamMixedModel, to get coefficiennt from RunVoomLimma use GetContrastResults
#' @param coefficient.number corresponds to the contrast, the nmber is in order of the contrast matrix
#' @param contrast.name this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.
#' @return a list of dataframes with contrast results indexed by cell type
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate full_join select
#' @importFrom stats p.adjust
#' @importFrom plyr mapvalues
#' @export
#'
#' @examples
#'\dontrun{
#'fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit,
#'                                         coefficient.number = 1,
#'                                         contrast.name = "foldchangedifference")
#'
#' }
GetContrastResultsRaw =  function(limma.fit.object.list, coefficient.number, contrast.name){
  print('this function will be deprecated by ExtractResult')
  print("this function returns results from dreamMixedModel, to get coefficient from RunVoomLimma, use GetContrastResults
        GetContrastResults uses emperican Bayes shrinkage see https://github.com/GabrielHoffman/variancePartition/issues/4 ")
  ## pt 1 return ONLY gene, logFC, AveExpr, contrast, celltype from eBayes call to format data and add raw p values from lme4, correct with BH.
  pvals = lapply(limma.fit.object.list, function(x){ data.frame(P.Value = x$p.value[ ,coefficient.number], gene = rownames(x$p.value))})
  lapply(pvals, function(x) rownames(x) = NULL)

  # parameters to run ordinary t statistic
  coef = lapply(limma.fit.object.list, function(x){ x$coefficients[ ,coefficient.number] })
  stdev = lapply(limma.fit.object.list, function(x){ x$stdev.unscaled[ ,coefficient.number] })
  sigma_ = lapply(limma.fit.object.list, function(x){ x$sigma })
  # Note eBayes t statistics are NOT used, only for formatting output, next section adds raw p values and logFC returned by lme4
  ebf = lapply(limma.fit.object.list, function(x) suppressWarnings(eBayes(x))) ##
  contrast_result = GetContrastResults(limma.fit.object.list = ebf,
                                       coefficient.number = coefficient.number,
                                       contrast.name = contrast.name)
  contrast_result = lapply(contrast_result, function(x){ x %>% dplyr::select(gene, logFC, AveExpr, contrast, celltype) })
  result = t_statistic = list()
  for (i in 1:length(limma.fit.object.list)) {
    # compute ordinary t-statistic (see e.g. limma eBayes documentation: ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
    t_statistic[[i]] = data.frame(t = coef[[i]] / stdev[[i]] / sigma_[[i]]) %>% tibble::rownames_to_column("gene")
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
#' @return a matrix
#' @importFrom dplyr filter select bind_rows
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#'\dontrun{
#' le = scglmmr::GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.1, NES.filter = -Inf)
#' genesub = do.call(rbind, le) %$% gene %>% unique
#' mtx2 = scglmmr::GetGeneMatrix(result.list = res,
#'                               stat_for_matrix = "logFC",
#'                               gene_subset = genesub,
#'                               pvalfilter = -Inf,
#'                               logfcfilter = 0.1)
#'
#' pheatmap::pheatmap(mtx2,
#'                    breaks =seq(from = 0, to = 2,length.out = 99),
#'                    filename = paste0(figpath,"LEgenes_heatmap.pdf"))
#' }
GetGeneMatrix = function(result.list, gene_subset = NULL, stat_for_matrix = "logFC", pvalfilter, logfcfilter){

  if (is.null(gene_subset)) {
    gene_subset = do.call(rbind, result.list)$gene %>% unique
  } else{
    gene_subset = gene_subset
  }
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
#' @param coefficient.number the coefficient from the custom contrasts , check with head(result@coefficients)
#' @param contrast.name this can be arbitrary and does not have to match the result coefficient name but is designed to force user to know which coefficient they are using from the fitted contrast model.
#'
#' @return list of gene ranks by t stat -- use as argument to `RunFgseaOnRankList`
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr arrange select
#' @importFrom limma topTable
#' @export
#'
#' @examples
#'\dontrun{
#'test = scglmmr::GetRankResults(limma.fit.object.list = bl, coefficient.number = 1, "test")
#' }
GetRankResults = function(limma.fit.object.list, coefficient.number, contrast.name){
  print('this function will be deprecated by ExtractResult')
  print(" returning ranks based on emperical bayes moderated t statistic")
  ranks = list()
  for (i in 1:length(limma.fit.object.list)) {
    test = as.data.frame(limma::topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf))
    ranks[[i]] =
      test %>%
      tibble::rownames_to_column("gene") %>%
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


#' GetRankResultsRaw get list of gene ranks by raw t statistic for fGSEA. for results returned by dreamMixedModel
#'
#' @param contrast.result.raw.list results returned by dreamMixedModel
#'
#' @return list of gene ranks by t stat -- use as argument to `RunFgseaOnRankList`
#' @importFrom dplyr arrange select desc
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#'\dontrun{
#' fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit,
#'                                          coefficient.number = 1,
#'                                          contrast.name = "foldchangedifference")
#' fit_rank = scglmmr::GetRankResultsRaw(contrast.result.raw.list = fit_res)
#' }
GetRankResultsRaw = function(contrast.result.raw.list){
  print('this function will be deprecated by ExtractResult')
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


#' FgseaList - wrapper around fast gene set enrichment analysis with the fgsea R package https://bioconductor.org/packages/release/bioc/html/fgsea.html to implement on a list of ranks indexec by cell type.
#'
#' @param rank.list.celltype results returned by GetRankResultsRaw or GetRankResults
#' @param pathways modules / gene sets as a named list each a single vector of unique gene IDS
#' @param maxSize see fgsea package
#' @param minSize see fgsea package
#'
#' @return results from fgsea package indexed by celltype
#' @importFrom fgsea fgsea
#' @importFrom dplyr arrange filter
#' @export
#'
#' @examples
#'\dontrun{
#' t1hvl_rank = GetRankResultsRaw(limma.fit.object.list  = dreamfit,
#' coefficient.number = 1,
#' contrast.name = "contrastName")
#' register(SnowParam(4))
#' pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
#' gsealist = FgseaList(rank.list.celltype = t1hvl_rank, pathways = btm,  BPPARAM = pparam)
#' }
#' # usage:
FgseaList = function(..., rank.list.celltype, pathways,
                     maxSize = 500, minSize = 9) {
  print(paste0(' fgsea: parallelize with `BPPARAM` '))
  # init storage
  gsea = list()
  ct = names(rank.list.celltype)

  # run fgsea on all pathways for each list element
  for (i in 1:length(rank.list.celltype)) {
    print(paste0('fgsea: ', ct[i], " ", i, ' of ', length(ct),
                 ' Testing ', length(pathways), ' pathways '))
    gsea[[i]] = fgsea::fgsea(...,
                             stats = rank.list.celltype[[i]],
                             pathways = pathways,
                             maxSize = maxSize,
                             minSize = minSize) %>%
      dplyr::mutate(celltype = ct[i]) %>%
      dplyr::arrange(pval)
  }
  names(gsea) = ct
  return(gsea)
}

#' LeadingEdgeIndexed - extract leading edge genes from FgseaList results, return a new embedded list of named modules, indexed by celltype
#' @param gsea.result.list results from FgseaList or RunFgseaOnRankList
#' @param padj.threshold within each cell type return list of leading edge genes with padj less than this parameter value.
#' @return embedded list - each list level 1 indexed by cell type, contains a new list level 2, the leading edge genes from the gsea results filtered by padj.
#' @importFrom dplyr filter select
#' @export
#'
#' @examples
#'\dontrun{
#' t1hvl_rank = GetRankResultsRaw(limma.fit.object.list  = ebf,
#' coefficient.number = 1,
#' contrast.name = "time_1_highvslow")
#' gsea = FgseaList(rank.list.celltype = t1hvl_rank)
#' celltype.indexed.modules.leadingedge = LeadingEdgeIndexed(gsea, 0.05)
#' }
LeadingEdgeIndexed = function(gsea.result.list, padj.threshold = 0.05){
  x = gsea.result.list
  newlist = list()
  for (i in 1:length(x)) {
    y = x[[i]] %>%
      dplyr::filter(padj < padj.threshold) %>%
      dplyr::select(pathway, leadingEdge)
    newlist[[i]] = y$leadingEdge
    names(newlist[[i]]) = y$pathway
  }
  names(newlist) = names(x)
  return(newlist)
}

#' EnrichmentJaccard - using gsea list and LeadingEdgeIndexed result, compute pairwise jaccard index of leadingedge genes within celltypes. saves a heatmap of modules for each cell type in savpath if saveplot = TRUE. Returns a gsea result dataframe with all celltypes combined and module annotated with average within celltype jaccard index and leadingedge genes.
#' @param gsealist results from FgseaList or RunFgseaOnRankList (recommend first lapply filter(padj < 0.05 e.g.) )
#' @param indexedgenes results fro mLeadingEdgeIndexed
#' @param saveplot if TRUE saves jaccard index heatmap to figpath
#' @param figpath place to save figures, a file.path().
#' @return curated dataframe of gsea results with average jaccard index.
#' @importFrom pheatmap pheatmap
#' @importFrom GeneOverlap newGOM getMatrix
#' @export
#'
#' @examples
#'\dontrun{
#'# read baseline enrichemnt results
#' g0 = FgseaList(rank.list.celltype = t1hvl_rank, pathways = btm,  BPPARAM = pparam)
#' filtered_g0 = lapply(g0, function(x) x %>% filter(padj < 0.05))
#'
#'   compute jaccard index of leadingedge genes within celltype
#'   li = LeadingEdgeIndexed(gsea.result.list = g0,padj.threshold = 0.05)
#'
#'   # enrichment jaccard
#'   d = EnrichmentJaccard(gsealist = filtered_g0, indexedgenes = li,
#'   saveplot = TRUE, figpath = figpath,
#'   fontsize_row = 7.5, fontsize_col = 7.5)
#'   d = d %>%
#'   mutate(leadingEdge = map_chr(leadingEdge, toString)) %>%
#'   select(celltype, av_jaccard,everything())
#'   write_delim(d,file = paste0(datapath, 'g0jaccard.csv'),delim = ',')
#' }
#'
EnrichmentJaccard = function(..., gsealist, indexedgenes, saveplot = FALSE, returnJaccardMtx = FALSE, figpath){

  # dat input check
  # remove enrichments from cell types without more than 2 modules enriched
  # add these back at end
  gsealist2 = gsealist
  subs = lapply(gsealist,nrow) > 2
  gsealist = gsealist[subs]
  indexedgenes = indexedgenes[subs]
  # confirm order of celltypes (lists) are the same
  stopifnot(isTRUE(all.equal(names(gsealist), names(indexedgenes))))
  jmat = list()
  for (i in 1:length(indexedgenes)) {
    print(names(indexedgenes)[i])
    #calculate pairwise jaccard index matrix
    overlap = GeneOverlap::newGOM(gsetA = indexedgenes[[i]], gsetB = indexedgenes[[i]], genome.size = NULL)
    jaccard_matrix = GeneOverlap::getMatrix(object = overlap, name = 'Jaccard')
    mean_ji = rowMeans(jaccard_matrix)
    if (isTRUE(saveplot)) {
      ph = pheatmap::pheatmap(jaccard_matrix, silent = TRUE, clustering_method = 'complete',
                              fontsize_col = 5, fontsize_row = 5, width = 15, height = 15,
                              filename = paste0(figpath, names(indexedgenes)[i], '.pdf'))
    } else {
      ph = pheatmap::pheatmap(jaccard_matrix, silent = TRUE, clustering_method = 'complete')
    }
    jmat[[i]] = ph
    #cluster modules
    clustered_mods = ph$tree_row$labels[ph$tree_row$order]
    gsealist[[i]] = gsealist[[i]][match(clustered_mods, gsealist[[i]]$pathway), ]
    gsealist[[i]]$av_jaccard = mean_ji
  }
  # return format
  d2 = do.call(rbind, gsealist2[names(subs[subs==FALSE])])
  d = do.call(rbind,gsealist)
  d = rbind(d,d2, fill=TRUE)
  # format jaccard matrix pheatmap output
  names(jmat) = names(gsealist)
  if(isTRUE(returnJaccardMtx)){
    ret = list('sortedgsea' = d, 'jaccard_matrix_list' = jmat)
    return(ret)
    } else{
  return(d)
    }
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

  .Deprecated("FgseaList")

  #require(fgsea)
  gsea = list()
  for (i in 1:length(rank.list.celltype)) {
    gsea[[i]] =
      suppressMessages(fgsea::fgsea(pathways = pathways,
                                    stats = rank.list.celltype[[i]],
                                    maxSize = maxSize,
                                    minSize = minSize,
                                    nperm = nperm)) %>%
      dplyr::mutate(celltype = names(rank.list.celltype)[i]) %>%
      dplyr::arrange(pval)
    if (positive.enrich.only == TRUE) {
      gsea[[i]] = gsea[[i]] %>% dplyr::filter(ES > 0)
    }
  }
  names(gsea) = names(rank.list.celltype)
  return(gsea)
}




#' RbindGseaResultList - prepare gsea result list for visualization funcitons GseaBubblePlot or GseaBarPlot; called by PlotFgseaList
#'
#' @param gsea_result_list result returned by RunFgseaOnRankList
#' @param NES_filter filter out results below this NES
#' @param padj_filter filter out results above this adjusted p threshold
#'
#' @return a dataframe of subsetted gsea results for all celltypes
#' @importFrom dplyr select filter mutate
#' @export
#'
#' @examples
#'\dontrun{
#'d = scglmmr::RbindGseaResultList(gsea_result_list = gsea1,NES_filter = -Inf,padj_filter = 0.2)
#' }
RbindGseaResultList = function(gsea_result_list, NES_filter = -Inf, padj_filter = 0.1){
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
#' @param rbind_gsea_result_dataframe result returned from RbindGseaResultList
#' @param celltype_name name of celltype to be plotted
#' @param save_path file path to save results
#' @param title title of plot
#' @param save_name name of file saved to save_path
#' @param fill_color color of bar
#' @param width ggsave param
#' @param height ggsave param
#'
#' @return nothing
#' @import ggplot2
#' @export
#'
#' @examples
#'\dontrun{
#'scglmmr::GSEABarPlot(d, celltype_name = 'celltype1', save_path = figpath, title = 'ct2', fill_color = 'dodgerblue' save_name = "plot.pdf")
#' }
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
#' @param save_path file path to save results
#' @param include_negative TRUE/FALSE whether to include negative enrichment in the plot.
#' @param save_name name of file saved to save_path
#' @param width ggpsave param
#' @param height ggsave param
#'
#' @return nothing
#' @import ggplot2
#' @export
#'
#' @examples
#'\dontrun{
#'scglmmr::GSEABubblePlot(d, save_path = figpath, save_name = "plot.pdf")
#' }
GSEABubblePlot = function(rbind_gsea_result_dataframe, save_path,  include_negative = TRUE, save_name, width = 8.5, height = 7.2) {

  # apply same aes for both includeneg with exception
  plot_param = list (
    geom_point(shape = 21),
    theme_bw(),
    scale_x_discrete(position = "top"),
    theme(axis.text.x=element_text(angle = 45, hjust = 0)),
    theme(axis.title.y = element_blank()),
    labs(fill = 'Normalized \n Enrichment \n Score', size = '-log10(padj)'),
    theme(legend.title = element_text(colour = "black", size = 8)),
    theme(axis.text.y = element_text(size = 8, color = "black")),
    theme(axis.text.x = element_text(size = 8.5, color = "black")),
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
  if (isTRUE(returnplot)) {
    return(p)
  }
  ggsave(p, filename = paste0(save_path,save_name,".pdf"), width = width, height = height)
  print(p)
  }
}


#' PlotFgsea identical to GSEABubblePlot, returns plot for manual adjustment or saving and also clusters the map.
#'
#' @param rbind_gsea_result_dataframe dataframe returned by RbindGseaResultList
#' @return ggplot object
#' @import ggplot2
#' @importFrom dplyr select
#' @importFrom tidyr spread
#' @importFrom tibble column_to_rownames
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @examples
#'\dontrun{
#'scglmmr::GSEABubblePlot(d, save_path = figpath, save_name = "plot.pdf")
#' }
PlotFgsea = function(gsea_result_list, NES_filter = -Inf, padj_filter = 0.1) {

  # combine into dataframe
  d = RbindGseaResultList(gsea_result_list, NES_filter = -Inf, padj_filter = 0.1)

  # cluster results
  hcdat = d %>%
    dplyr::select(celltype, pathway,NES) %>%
    tidyr::spread(celltype, NES) %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  hcdat[is.na(hcdat)] = 0
  xx = pheatmap::pheatmap(hcdat, silent = TRUE, clustering_method = "average")
  module_order = xx$tree_row$labels[xx$tree_row$order]
  celltype_order = xx$tree_col$labels[xx$tree_col$order]
  d$celltype = factor(d$celltype, levels = celltype_order)
  d$pathway = factor(d$pathway, levels = module_order)


  # plot aes
  plot_param = list (
    geom_point(shape = 21),
    theme_bw(),
    scale_x_discrete(position = "top"),
    theme(axis.text.x=element_text(angle = 45, hjust = 0)),
    theme(axis.title.y = element_blank()),
    labs(fill = 'Normalized \n Enrichment \n Score', size = '-log10(padj)'),
    theme(legend.title = element_text(colour = "black", size = 8)),
    theme(axis.text.y = element_text(size = 8, color = "black")),
    theme(axis.text.x = element_text(size = 8.5, color = "black")),
    guides(shape = guide_legend(override.aes = list(size = 5))),
    guides(color = guide_legend(override.aes = list(size = 5))),
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  )
  p = ggplot(d, aes(y = pathway, x = celltype, fill = NES, size = n_logp)) +
    plot_param + xlab("") +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0)
    return(p)
}


#' GetLeadingEdgeFull Get a tidy dataframe of ALL Leading Edge Genes from gene set enrichment for all cell types
#'
#' @param gsea.list results returned by RunFgseaOnRankList
#' @param padj.filter filter reslts
#' @param NES.filter filter results
#'
#' @return a list
#' @importFrom purrr map_int
#' @importFrom dplyr filter select
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#'\dontrun{
#'lefull = scglmmr::GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.1,NES.filter = -Inf)
#' }

GetLeadingEdgeFull = function(gsea.list, padj.filter, NES.filter){
  gseasub = lapply(gsea.list, function(x){
    x %>%
    dplyr::filter(NES > NES.filter & padj < padj.filter) %>%
    dplyr::select(pathway, leadingEdge)
  })
  # get lists with results passing filter and subset result lists by that index
  g =  which(purrr::map_int(gseasub, nrow) > 0) %>% names
  gseasub = gseasub[g]
  stopifnot(g > 0 ); print("lists passing filter") ; print(g)
  mods = lapply(gseasub, function(u){
    testcase = u; paths = u$pathway
    dataret = data.frame()

    for (i in 1:length(paths)) {
      genes_use = testcase %>% dplyr::filter(pathway == paths[i])
      genes_use = genes_use$leadingEdge %>%
        unlist() %>%
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
#' @return a vector of genes
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#'\dontrun{
#' le_mono = GetLeadingEdgeGenes(gsea.result.list = gsea1, celltype.index = 4, module.name = 'my_modulename_from_gsearesults')
#' }
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


#' CombineResults - For all cell types merge the gsea leading edge genes with their contrast model coefficieint and p value from limma / dream
#'
#' @param gsealist list of results returned by `RunFgseaOnRankList()`
#' @param contrastlist list of results returned by `GetContrastResults()` or `GetContrastResultsRaw()`
#' @param gseafdr the adjusted p value threshold to filter out results (gsea)
#' @param genefdr the adjusted p value threshold to filter out individual genes recommend keep this high.
#'
#' @return a tidy dataframe
#' @importFrom dplyr mutate filter select group_by
#' @export
#'
#' @examples
#'\dontrun{
#' combined_results = CombineResults(gsealist = testgsea, contrastlist = testmod, gseafdr = 0.05,genefdr = 0.2)
#' }
CombineResults = function(gsealist, contrastlist, gseafdr, genefdr){


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

  # remove any celltypes without gsea results passing filter
  enriched_dim = sapply(gs, dim)
  remove = which(enriched_dim[1, ] == 0) %>% names

  # subset the DE gene list by the celltypes with BTM results
  generes = generes %>% dplyr::filter(!celltype %in% remove)

  # subset the gsea list by the celltypes with DE genes.
  celltype_result = celltype_result[!celltype_result %in% remove]
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
#' @param av.exprs.list object returned by `PseudobulkList` summed or average counts
#' @param gsea.list object returned by RunFgseaOnRankList
#' @param padj.filter filter for adjusted p from GSEA
#' @param NES.filter filter for normalized enrichment score from GSEA
#'
#' @return a list of tidy dataframes by celltype
#' @importFrom purrr map_int
#' @importFrom dplyr filter select
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather everything
#' @export
#'
#' @examples
#'\dontrun{
# make tidy average data for visualization of weighted pb results
#' av = scglmmr::PseudobulkList(rawcounts = umi,
#'                              metadata = meta,
#'                              sample_col = "sample",
#'                              celltype_col = "celltype",
#'                              avg_or_sum = 'average')
#' le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av,
#'                                            gsea.list = hlmk_ctm0,
#'                                            padj.filter = 0.1,
#'                                            NES.filter = -Inf)
#' }
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
      base::as.data.frame() %>%
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
#' @param av.exprs.list object returned by `PseudobulkList` summed or average counts
#' @param result.list object returned by object returned by `GetContrastResultsRaw()` or `GetContrastResults()`
#' @param P.Value.filter filter results
#' @param logFC.filter filter results
#' @param top_n_genes instead of stat thresholds above, get expression of n genes specified by this param, subsets each celltype by the top n genes  ranked by t statistic.
#'
#' @return a list of tidy dataframes by celltype
#' @importFrom dplyr filter arrange
#' @importFrom purrr map_int
#' @importFrom tibble rownames_to_column
#' @tidyr gather
#' @export
#'
#' @examples
#'\dontrun{
# make tidy average data for visualization of weighted pb results
#' av = scglmmr::PseudobulkList(rawcounts = umi,
#'                              metadata = meta,
#'                              sample_col = "sample",
#'                              celltype_col = "celltype",
#'                              avg_or_sum = 'average')
#' fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit,
#'                                         coefficient.number = 1,
#'                                         contrast.name = "foldchangedifference")
#' le_expr = scglmmr::TopGenesTidySampleExprs(av.exprs.list = av, result.list = fit_res, P.Value.filter = 0.2, logFC.filter=0.1, top_n_genes = 20)
#'}
TopGenesTidySampleExprs = function(av.exprs.list, result.list, P.Value.filter, logFC.filter, top_n_genes = NULL){

  resultsub = lapply(result.list, function(x){
    x = x %>%
      dplyr::filter(logFC > logFC.filter & P.Value < P.Value.filter) %>%
      dplyr::arrange(logFC) %$%
      gene
  })

  # # get top n genes if specified.
  if (!is.null(top_n_genes)){
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
      suppressMessages(tidyr::gather(sample, av_exp, index1:index2))

    tidy[[i]]$celltype = names(avsub[i])
  }
  names(tidy) = names(avsub)
  return(tidy)
}

#' LeadEdgeSampleHeatmap make a heatmap of average expression for top or leading edge genes returned by LeadEdgeTidySampleExprs
#'
#' @param tidy.exprs.list the object returned by returned by TopGenesTidySampleExprs or LeadEdgeTidySampleExprs which is created from metadata adn the object returned by PseudobulkList (can be average or summed expression)
#' @param modulename The name of the module to plot
#' @param celltype_plot the name of the celltype to plot
#' @param metadata cells x metadata dataframe
#' @param metadata_annotate a vector of variables (columns) of sample level  metadata to annotate on the heatmap. Must be categorical for each sample e.g. "age" but not "nUMI"
#' @param sample_column the column in the metadata object corresponding to the sample labels, usually 'smaple'
#' @param returnmat instead of making an annotated heatmap just return the matrix of averge values per sample of the module subset
#'
#' @return a pheatmap object
#' @importFrom dplyr filter select group_by summarize_each
#' @importFrom tidyr spread
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @examples
#'\dontrun{
#' scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr,
#'                                modulename = "MODULENAME",
#'                                elltype_plot = "TCELL",
#'                                metadata = meta,
#'                                metadata_annotate = c('group', 'timepoint', 'age', 'gender'),
#'                                sample_column = 'sample',
#'                               returnmat = F,
#'                                savepath = figpath,
#'                               savename = "filename")
#'}
LeadEdgeSampleHeatmap = function(tidy.exprs.list, modulename, celltype_plot,
                                 metadata, metadata_annotate, sample_column, returnmat = FALSE,
                                 plotwidth = 5, plotheight = 8, savepath , savename ){

  ### subset average expression object
  d = tidy.exprs.list[[celltype_plot]] %>%
    dplyr::filter(module == modulename) %>%
    dplyr::select(gene, sample, av_exp) %>%
    tidyr::spread(sample, av_exp) %>%
    tibble::column_to_rownames("gene")
  if(isTRUE(returnmat)) {
    return(d)
  } else{
    gvar = rlang::sym(sample_column)
    heatmap_anno = meta[meta$celltype == celltype_plot,   c(sample_column, metadata_annotate)] %>%
      dplyr::group_by({{gvar}}) %>%
      dplyr::summarise_each(list(~unique(.))) %>%
      tibble::column_to_rownames(sample_column)

    # cu = rev(pals::brewer.rdbu(12))
    cu = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3",
           "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")
    x = pheatmap::pheatmap(d, scale = 'row',
                           annotation = heatmap_anno,
                           color = cu,
                           width = plotwidth, height = plotheight,
                           filename = paste0(savepath, savename, ".pdf"))
    return(x)
  }
}



#' RunHypergeometricTest - run a hypergeometric test on results returned by GetContrastResults or GetContrastResultsRaw
#'
#' @param result_list results returned by by GetContrastResults or GetContrastResultsRaw
#' @param TERM2GENE_dataframe see clusterprofiler, these objects are automatically loaded in the scglmmr package one of term_df_btm, term_df_kegg, term_df_reactome etc.
#' @param pval_threshold p threshold for genes to consider in hypergeometric distribution
#' @param logFC_threshold logFC threshold  for genes to consider in hypergeometric distribution
#' @param usefdr_threshold use the FDR adjusted p values for ranking genes-this is a strict filter for single cell data, recommended to set FALSE
#'
#' @return tidy hypergeometric test results dataframe
#' @importFrom AnnotationDbi keys
#' @importFrom clusterProfiler enricher
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr mutate bind_rows filter
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#'\dontrun{
### Hypergeometric erichment
#' load(termdf) # this term2gene dataframe is included in the package see clusterProfiler
#' hyp = scglmmr::RunHypergeometricTest(result_list = fit_res,
#'                                      TERM2GENE_dataframe = termdf,
#'                                      pval_threshold = 0.1,
#'                                      logFC_threshold = 0,
#'                                      usefdr_threshold = FALSE)
#' # plot results
#' scglmmr::PlotHypergeometric(hyperg_result = hyp,
#'                             p.adjust.filter = 0.1,
#'                             genenumber_filter = 2,
#'                             savepath = figpath,
#'                             savename = "name",
#'                             title = "title")
#'  }
RunHypergeometricTest = function(result_list, TERM2GENE_dataframe, pval_threshold = 0.05,
                                 logFC_threshold = 0.5, usefdr_threshold = FALSE){

  print("result_list = list of dataframes indexed by celltype \n column names: gene, logFC, adj.P.Val, P.Value; use PlotHypergeometric() on results")

  # create list of entrez IDs for enriched features passing specified filters
  entrez_subset = list()
  for (i in 1:length(result_list)) {
    if (isTRUE(usefdr_threshold)) {
      d = result_list[[i]]
      d = d[d$adj.P.Val < pval_threshold & d$logFC > logFC_threshold, ]
    } else {
      d = result_list[[i]]
      d = d[d$p_val < pval_threshold & d$logFC > logFC_threshold, ]
    }
    if (nrow(d) < 1) {
      entrez_subset[[i]] = d = NA
    } else {
      # map geneID to entrez ids
      ent =
        tryCatch(
          AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                keys = d$gene, columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL"),
          error = function(e) return(NA)
        )
      # if after mapping, no ids are mapped, keep the NA value.
      if (all(is.na(ent))) {
        entrez_subset[[i]] = NA
      } else {
        dup = duplicated(ent$SYMBOL)
        ent = ent[!dup, ]$ENTREZID
        ent = ent[!is.na(ent)]
        entrez_subset[[i]] = ent
      }
    }
  }
  names(entrez_subset) = names(result_list)
  # remove individual subsets without any unmapped features
  entrez_subset = entrez_subset[!is.na(entrez_subset)]
  # run hypergeometric test
  # init strage and iterate over remaining subsets
  hypergeometric = list()
  for (i in 1:length(entrez_subset)){

    cells = names(entrez_subset[i])
    print(paste0("hypergeometric test in: ", cells, " index ", i))
    # don't run hypergeometric test if there is only 1 gene enriched or 0 genes enriched.
    if(length(entrez_subset[[i]]) <= 1){
      hypergeometric[[i]] = NA
    } else {
      # run hypergeometric test with clusterprfiler package
      hypergeometric[[i]] =
        suppressMessages(
        tryCatch(
          clusterProfiler::enricher(entrez_subset[[i]], TERM2GENE = TERM2GENE_dataframe)@result %>%
            dplyr::mutate(celltype = cells) %>%
            tidyr::separate(GeneRatio,into = c("gene_num", "gene_denom"), sep = "/") %>%
            dplyr::mutate(gene_num = as.numeric(gene_num)) %>%
            dplyr::mutate(gene_denom = as.numeric(gene_denom)) %>%
            dplyr::mutate(gene_ratio = round(gene_num / gene_denom, 3)),
          error = function(e) return(NA)
        )
        )
    }
  }
  # combine results and format for PlotHypergeometric() function
  hypergeometric = hypergeometric[!is.na(hypergeometric)]
  hyp = hypergeometric %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(n_logp = -log10(p.adjust))
  return(hyp)
}



#' PlotHypergeometric - plot results returned by RunHypergeometricTest
#'
#' @param hyperg_result result returned by `RunHypergeometricTest`
#' @param p.adjust.filter filter results
#' @param genenumber_filter number of genes within an enrichment minimum.
#' @param savepath save oath
#' @param savename name of object saved to save_path
#' @param title title of plot
#' @param height ggsave param
#' @param width ggsave param
#'
#' @return nothing
#' @import ggplot2
#' @export
#'
#' @examples
#'\dontrun{
### Hypergeometric erichment
#' load(termdf) # this term2gene dataframe is included in the package see clusterProfiler
#' hyp = scglmmr::RunHypergeometricTest(result_list = fit_res,
#'                                      TERM2GENE_dataframe = termdf,
#'                                      pval_threshold = 0.1,
#'                                      logFC_threshold = 0,
#'                                      usefdr_threshold = FALSE)
#' # plot results
#' scglmmr::PlotHypergeometric(hyperg_result = hyp,
#'                             p.adjust.filter = 0.1,
#'                             genenumber_filter = 2,
#'                             savepath = figpath,
#'                             savename = "name",
#'                             title = "title")
#'  }
PlotHypergeometric = function(hyperg_result, p.adjust.filter = 0.05, genenumber_filter = 0,
                              savepath = figpath, savename , title, height = 10, width = 8 ){

  # filter results
  hyp = hyperg_result[hyperg_result$gene_num > genenumber_filter & hyperg_result$p.adjust < p.adjust.filter, ]
  # plot
  p = ggplot(hyp, aes(y = ID, x = celltype, fill = n_logp, size = gene_ratio)) +
    geom_point(shape = 21) +
    scale_fill_viridis_c(option = "B") +
    theme_bw() +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
    theme(axis.title.y = element_blank()) +
    labs(fill = '-log10(ajdusted P value)', size = 'gene ratio') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) +
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) +
    guides(shape = guide_legend(override.aes = list(size = 5))) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    theme(plot.title = element_text(face = "bold", color = "black", size = 8)) +
    ggtitle(title)
  # save
  ggsave(p, filename = paste0(savepath, savename, ".pdf"), width = width, height = height)
}


##### plot average gene distributions for each sample in each cohort.
# most general version

#' GetTidySummary - tidy data summary for a single cell type
#'
#' @param av.exprs.list - object returned by PseudobulkList
#' @param celltype.index - index of celltype to ret results see names of PseudobulkList object
#' @param genes.use - subet of genes to use
#'
#' @return a tidy dataframe
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @export
#'
#' @examples
#' #' @examples
#'\dontrun{
#' gene_highlight =  c("IRF1","TNFRSF17","ABL1")
#' mono = GetTidyCohort(av.exprs.list = av, celltype.index = 7, genes.use = gene_highlight)
#' PlotGeneDistCohort(merged_av_data = mono,
#'                    save_name = "mono_highlight_512",
#'                    save_path = figpath,
#'                    title = paste0(names(av[7]), "genesub" ),
#'                    height = 3.8, width = 4.5,
#'                    nrow = 2)
#'  }
GetTidySummary = function(av.exprs.list, celltype.index, genes.use){
  gene.index.1 = genes.use[1]
  gene.index.2 = genes.use[length(genes.use)]
  tidy_data =
    av.exprs.list[[celltype.index]][genes.use, ] %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather(key = gene, value = count, gene.index.1:gene.index.2)
  return(tidy_data)
}



# plot the tidy gene x sample summary by cohort.
#' PlotGeneDistCohort
#'
#' @param merged_av_data data returned for a single cell type by `GetTidySummary`
#' @param save_path file path to save results
#' @param save_name name of plot saved to `save_path`
#' @param title title of plot
#' @param nrow number of rows to facet genes plotted on
#' @param height ggsave param
#' @param width ggsave param
#' @param plot_subset whether to subset to some of the genes in `merged_av_data`
#' @param genes_plot the subset of genes to plot
#'
#'@import ggplot2
#'
#' @return nothing
#' @export
#'
#' @examples
#' @examples
#'\dontrun{
#' gene_highlight =  c("IRF1","TNFRSF17","ABL1")
#' mono = GetTidyCohort(av.exprs.list = av, celltype.index = 7, genes.use = gene_highlight)
#' PlotGeneDistCohort(merged_av_data = mono,
#'                    save_name = "mono_highlight_512",
#'                    save_path = figpath,
#'                    title = paste0(names(av[7]), "genesub" ),
#'                    height = 3.8, width = 4.5,
#'                    nrow = 2)
#'  }
PlotGeneDistCohort = function(merged_av_data,
                              save_path, save_name,
                              title = NULL,
                              nrow = 5, height = 8, width = 10,
                              plot_subset = FALSE, genes_plot = NULL){

  cu = c("dodgerblue", "midnightblue", "red", "firebrick",  "#FFBB78FF", "#FF7F0EFF")
  if(plot_subset == TRUE) {
    p = ggplot(merged_av_data %>% filter(gene %in% genes_plot), aes(x = group, y = count, fill = interaction(timepoint, group ))) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free_y", nrow = nrow) +
      theme_bw(base_size = 10.5)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = cu) +
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
      scale_fill_manual(values = cu) +
      theme(strip.background = element_blank()) +
      theme(strip.text = element_text(face = "bold",family = "Helvetica")) +
      theme(axis.text.x =  element_blank()) +
      theme(axis.text.y =  element_text(size = 6)) +
      ggtitle(title)+
      ylab("Average Expression")
    ggsave(p, filename = paste0(save_path,save_name,".pdf"), height = height,  width = width)
  }
}


# retired util functions

# ### match order of 2 lists of data frames based on "celltype" column useful to e.g. match separately run models by fgsea result lists
# MatchfgseaResultIndex = function(list_to_reference, list_to_reorder){
#
#   print("cell type names must match")
#
#   result.list1 = list_to_reference
#   result.list2 = list_to_reorder
#
#   ct = sapply(result.list1, function(x){ unique(x$celltype)}) %>% unname
#   print("list 1 order of results by cell type")
#   print(ct)
#   ct2 = sapply(result.list2, function(x){ unique(x$celltype) }) %>% unname
#   print("list 2 order of results by cell type")
#   print(ct2)
#
#   if (!(length(ct) == length(ct2))) {
#     print("the celltype result lists are unequal length and will be merged by the intersection")
#   }
#
#   ct2 = ct2[ct2 %in% ct]
#   ct = ct[ct %in% ct2]
#
#   reorder_2_match1 = ct2[order(match(ct2, ct))]
#
#   result.list2 = result.list2[reorder_2_match1]
#   return(result.list2)
# }

# # Make volcano plot for each celltype
# VolcanoPlotTop = function(contrast.result.list, contrast.name, save.path, size = 3, fig.height, fig.width) {
#   require(ggrepel)
#   for (i in 1:length(contrast.result.list)) {
#     #format result list for volcano plot
#     contrast.result.list[[i]] = contrast.result.list[[i]] %>%
#       arrange(adj.P.Val)  %>%
#       mutate(nlp = -log10(adj.P.Val)) %>%
#       mutate(color =
#       if_else(nlp > 1 & LogFC > 0.2, true = "1",
#       if_else(nlp > 1 & LogFC < -0.2,  true = "2",
#       if_else(nlp > 1 & LogFC < 0.2 & LogFC > - 0.2,  true = "3", false = "0"))))
#
#     # volcano plot
#     p = ggplot(data = contrast.result.list[[i]], aes(x = LogFC, y = nlp, label = gene, color = color)) +
#       geom_point(show.legend = FALSE) +
#       theme_bw() +
#       geom_text_repel(data= contrast.result.list[[i]] %>%
#                         filter(nlp > 1.3) %>%
#                         filter(abs(LogFC) > 0.2),
#                       aes(label=gene, color = color), size = size ,show.legend = FALSE) +
#       geom_vline(xintercept = 0.2, linetype = "dashed") +
#       geom_vline(xintercept =  - 0.2, linetype = "dashed") +
#       geom_hline(yintercept = 1.3, linetype = "dashed") +
#       labs(x = "Log FC", y = "-log10 adjusted p value") +
#       scale_color_manual(values = c("black", "red", "blue", "black")) +
#       theme(text =  element_text(colour = "black", family = "Helvetica",size = 12)) +
#       ggtitle(names(contrast.result.list[i]))
#
#     # save plot
#     ggsave(plot = p,
#            filename = paste0(contrast.name, names(contrast.result.list)[i],".pdf"),
#            path = save.path,
#            width = fig.width, height = fig.height)
#   }
# }
#
#
