# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com


#' ExtractResult - convenience function to return statistics for downstream analysis functions such as FgseaList. Returns results from list of dream or lmFit results from those functions natively (use model.fit.list = list(fit)) or from scglmmr::RunVoomLimma and scglmmr::dreamMixedModel
#' @param model.fit.list list of model results indexed by celltypes returned by `scglmmr::dreamMixedModel`,  `scglmmr::RunVoomLimma`, or manually by `lmFit` or `dream`.
#' @param what what to return what = c('statistics', 'lmer.z.ranks' or 'gene.t.ranks')[1] defaults to statistics for each cell type, e.g. avg exprs, logFC, t statistic, pval, adj.P.Val etc. If gene.z.ranks, ranks genes based on z statistic (mixed models) and returns a list (indexed by celltype of named numeric vector of genes ranked for FgseaList.
#' @param coefficient.number what coefficient to return -- this needs to be one of model.fit.list$coefficients: check the order of the coefficients. Results returned from dream include statistical contrasts and estimated coefficients from the model. If limma::contrasts.fit was used (e.g. if using do_contrast_fit = TRUE in RunVoomLimma), these 'coefficients' are results of the statistical contrast.
#' @param coef.name the name of the estimated coefficient for which results are being returned; if returning results from a statistical contrast e.g. limma::contrasts.fit() this will be the name of the contrast. If returning a model fit with the dream function, can also be contrast specified by from variancePartition::makeContrastsDream() or a fixed effect parameter that was included in the model.
#' @return a list of dataframes with contrast results indexed by cell type or a list of genes ranked b t statistic in format ready for FgseaList.
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr mutate select
#' @export
ExtractResult = function(model.fit.list,
                         what = c('statistics', 'lmer.z.ranks', 'gene.t.ranks')[1],
                         coefficient.number,
                         coef.name) {
  values.possible = c('statistics', 'lmer.z.ranks', 'gene.t.ranks')
  if (!what %in% values.possible) {
    stop('`what` must be one of: statistics, lmer.z.ranks, gene.t.ranks')
  }


  #init
  celltypes = names(model.fit.list)

  # check user input data
  model_type = model.fit.list[[1]]$method
  message1 = 'lmer mixed model fit'
  message2 = 'linear model fixed effects only fit'
  messageprint = ifelse(model_type == 'lmer', message1, message2)
  print(paste0('returning results of ', model_type, ' model: ', messageprint))

  # ensure user estimating the correct coefficient
  coefs = colnames(model.fit.list[[1]]$coefficients)
  print('coefficients available from model fit object: ')
  print(coefs)
  stopifnot(all.equal(as.character(coef.name),
                      as.character(coefs[coefficient.number])))

  # print data returning
  if (what == 'statistics') {
    print('returning statistics for: ')
    print(coefs[coefficient.number])
    print('to return ranks change argument to `what`')
  } else{
    print('returning gene ranks for: ')
    print(coefs[coefficient.number])
  }

  # extract results
  test = list()
  if (model_type == 'lmer') {
    # use variancePartition topTable function for mixed model
    for (i in 1:length(model.fit.list)) {
      test[[i]] =
        variancePartition::topTable(fit = model.fit.list[[i]],
                                    coef = coefficient.number,
                                    number = Inf) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(contrast = rep(coef.name)) %>%
        dplyr::mutate(celltype = celltypes[i]) %>%
        arrange(desc(z.std)) # arrange by signed z statistic calculated by dream
    }
  } else {
    # use limma topTable function for fixed effects model
    for (i in 1:length(model.fit.list)) {
      test[[i]] =
        limma::topTable(fit =  model.fit.list[[i]],
                        coef = coefficient.number,
                        number = Inf) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(contrast = rep(coef.name)) %>%
        dplyr::mutate(celltype = celltypes[i]) %>%
        arrange(desc(t)) # arrange by empirical bayes moderated t stat
    }
  }
  # convenince for returning genes ranked (as above from the arrange call)
  # added the base::sort, decreasing = TRUE line to futureproof in case arrange(desc(z.std)) removed
  if (what == 'gene.t.ranks') {
    ret = lapply(test, function(x) {
      structure(x$t, names = as.character(x$gene))
    })
    ret = lapply(ret, sort, decreasing = TRUE)
  } else if (what == 'lmer.z.ranks') {
    ret = lapply(test, function(x) {
      structure(x$z.std, names = as.character(x$gene))
    })
    ret = lapply(ret, sort, decreasing = TRUE)
  } else{
    # default - return full model statistics
    ret = test
  }
  names(ret) = celltypes
  return(ret)
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


#' HeatmapDiag utility function for `pheatmap` using `slanter` to order the maximum values across the matrix diagonal
#'
#' @param matrix an R matrix as input to `pheatmap`. Could be the object returned by `scglmmr::GetGeneMatrix`
#' @return a pheatmap object
#' @importFrom slanter slanted_orders
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @examples
#'\dontrun{
#' # make gene plot of top 50 ranked genes within each subst
#' # fit1 = dreamMixedModel(...)
#' r1 = ExtractResult(model.fit.list = fit1, coefficient.number = 1,coef.name = 'L1', what = 'gene.t.ranks')
#' rank = lapply(r1,function(x) names(x[1:50]))
#' top50ranks = unlist(rank) %>% unique()
#' res1 = ExtractResult(model.fit.list = fit1, coefficient.number = 1,coef.name = 'L1')
#' mtx = GetGeneMatrix(result.list = res1, gene_subset = top50ranks, stat_for_matrix = 'logFC', pvalfilter = 1, logfcfilter = -Inf)
#' pdf(file = paste0(figpath, 'geneplot.pdf'), width = 3.5 ,height = 5)
#' HeatmapDiag(matrix = mtx, fontsize_row = 1)
#' dev.off()
#' }
HeatmapDiag = function(matrix, ...){
  mtx2 = matrix
  # rm neg only to calc diagonalizing order
  mtx2[mtx2<0] = 0
  orders = slanter::slanted_orders(data = mtx2, order_rows = TRUE, order_cols = TRUE)
  #plot the original data in diagonalized max intensity order
  roworder = rownames(matrix)[orders$rows]
  colorder = colnames(matrix)[orders$cols]
  p = pheatmap::pheatmap(matrix[roworder, colorder],cluster_rows = FALSE, cluster_cols = FALSE, ...)
  return(p)
}



#' TidySampleData convert data from PseudobulkList into a dataframe for each sample across cell types of the top differentially expressed genes for a contrast
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
TidySampleData = function(av.exprs.list, result.list, P.Value.filter, logFC.filter, top_n_genes = NULL){

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


##### plot average gene distributions for each sample in each cohort.
# most general version

#' GetTidySummary - tidy data summary for a single cell type
#'
#' @param av.exprs.list - object returned by PseudobulkList (use average or first convert summed counts to cpm)
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







