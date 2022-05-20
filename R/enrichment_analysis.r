# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

#' FgseaList - wrapper around fast gene set enrichment analysis with the fgsea R package https://bioconductor.org/packages/release/bioc/html/fgsea.html to implement on a list of ranks indexec by cell type.
#'
#' @param rank.list.celltype results returned by GetRankResultsRaw or GetRankResults
#' @param pathways modules / gene sets as a named list each a single vector of unique gene IDS
#' @param maxSize see fgsea package
#' @param minSize see fgsea package
#'
#' @return results from fgsea package indexed by celltype
#' @importFrom fgsea fgsea
#' @importFrom dplyr arrange mutate
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
  print(' If using this method, cite Korotkevich et. al. Biorxiv (2021) doi.org/10.1101/060012')
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
LeadingEdgeIndexed = function(gsea.result.list, padj.threshold = 0.05, p.threshold = NULL){
  x = gsea.result.list
  newlist = list()
  # enable p thresholding
  if (!is.null(p.threshold)) {
    for (i in 1:length(x)) {
      y = x[[i]] %>%
        dplyr::filter(pval < p.threshold) %>%
        dplyr::select(pathway, leadingEdge)
      newlist[[i]] = y$leadingEdge
      names(newlist[[i]]) = y$pathway
    }
  } else{
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
  # remove enrichments from cell types without more than 1 modules enriched
  # add these back at end
  gsealist2 = gsealist
  subs = lapply(gsealist,nrow) > 1
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
    jmat[[i]] = jaccard_matrix
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
    ret = list('sortedgsea' = d, 'jaccard_matrix_list' = jmat )
    return(ret)
  } else{
    return(d)
  }
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
RbindGseaResultList = function(gsea_result_list, NES_filter = -Inf, padj_filter = NULL, pval_filter = NULL){

  score = lapply(gsea_result_list, function(x) {
    x = x %>%
      dplyr::select(pathway, pval, padj, NES, celltype) %>%
      dplyr::filter(NES > NES_filter)
  })
  # combine and filter
  if (!is.null(pval_filter)) {
    score = do.call(rbind, score) %>%
      dplyr::filter(pval < pval_filter) %>%
      dplyr::mutate(n_logp = -log10(padj))
  } else{
    score = do.call(rbind, score) %>%
      dplyr::filter(padj < padj_filter) %>%
      dplyr::mutate(n_logp = -log10(padj))
  }
  return(score)
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
#' p = PlotFgsea(gsea_result_list = g1, padj_filter = 0.1)
#' }
PlotFgsea = function(gsea_result_list, NES_filter = -Inf, padj_filter = 0.1, p.threshold = NULL) {
  # combine results
  if (!is.null(p.threshold)) {
    d = RbindGseaResultList(gsea_result_list,
                            NES_filter = -Inf,
                            pval_filter = 0.05)
    p.legend = '-log10(p)'
  } else {
    # combine into dataframe
    d = RbindGseaResultList(gsea_result_list,
                            NES_filter = -Inf,
                            padj_filter = 0.1)
    p.legend = '-log10(padj)'
  }


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
  # reorder as factors for ggplot
  d$celltype = factor(d$celltype, levels = celltype_order)
  d$pathway = factor(d$pathway, levels = module_order)


  # plot aes
  plot_param = list (
    geom_point(shape = 21),
    theme_bw(),
    scale_x_discrete(position = "top"),
    theme(axis.text.x=element_text(angle = 45, hjust = 0)),
    theme(axis.title.y = element_blank()),
    labs(fill = 'Normalized \n Enrichment \n Score', size = p.legend),
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
#' @param p.filter raw p value filter -- padj.filter must not be specified
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

GetLeadingEdgeFull = function(gsea.list, padj.filter = NULL, NES.filter, p.filter= NULL){

  # filter results
  if (!is.null(p.filter)) {
    gseasub = lapply(gsea.list, function(x) {
      x %>%
        dplyr::filter(NES > NES.filter & pval < p.filter) %>%
        dplyr::select(pathway, leadingEdge)
    })

  } else {
    gseasub = lapply(gsea.list, function(x) {
      x %>%
        dplyr::filter(NES > NES.filter & padj < padj.filter) %>%
        dplyr::select(pathway, leadingEdge)
    })
  }
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
#' @param contrastlist list of results returned by  `scglmmr::ExtractResult()` or from older versions of scglmmr, the equivalent results from: `GetContrastResults()` , `GetContrastResultsRaw()`
#' @param gseafdr the adjusted p value threshold to filter out results (gsea)
#' @param gseap the p value threshold to filter out results (gsea); gseafdr must not be specified to use this option
#' @param genefdr the adjusted p value threshold to filter out individual genes - by default all the leading edge genes are returned (recommended)
#'
#' @return a tidy dataframe
#' @importFrom dplyr mutate filter select group_by
#' @export
#'
#' @examples
#'\dontrun{
#' combined_results = CombineResults(gsealist = testgsea, contrastlist = testmod, gseafdr = 0.05,genefdr = 0.2)
#' }
CombineResults = function(gsealist, contrastlist, gseafdr, gseap, genefdr = Inf){

  # filter results
  if (!is.null(gseap)) {
    # filter gsea esults and cat pathway and celltype
    gs = lapply(gsealist, function(x) {
      x = x %>%
        dplyr::filter(pval < gseap) %>%
        dplyr::mutate(name = paste(pathway, celltype, sep = "__"))
    })
  } else{
    gs = lapply(gsealist, function(x) {
      x = x %>%
        dplyr::filter(padj < gseafdr) %>%
        dplyr::mutate(name = paste(pathway, celltype, sep = "__"))
    })
  }

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
#' @param p.filter filter for p from GSEA padj.filter must not be specified to use this option.
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
LeadEdgeTidySampleExprs = function(av.exprs.list, gsea.list, padj.filter, p.filter = NULL, NES.filter){


  # filter results
  if (!is.null(gseap)) {
    gseasub = lapply(gsea.list, function(x){x = x %>%
      dplyr::filter(NES > NES.filter & pval < p.filter) %>%
      dplyr::select(pathway, leadingEdge)
    })
  } else{
    gseasub = lapply(gsea.list, function(x){x = x %>%
      dplyr::filter(NES > NES.filter & padj < padj.filter) %>%
      dplyr::select(pathway, leadingEdge)
    })
  }
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

#' LeadEdgeSampleHeatmap make a heatmap of average expression for top or leading edge genes returned by LeadEdgeTidySampleExprs
#'
#' @param tidy.exprs.list the object returned by returned by TidySampleData or LeadEdgeTidySampleExprs which is created from metadata and the object returned by PseudobulkList (average expression prefered for visualization; else recommend e.g. first lapply(x, edgeR::cpm) to standardize the bulk data.
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
#'                                celltype_plot = "TCELL",
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


#' GSEABarPlot - plot gsea results for a single cell type
#'
#' @param rbind_gsea_result_dataframe result returned from RbindGseaResultList
#' @param celltype_name name of celltype to be plotted
#' @param fill_color color of bar
#' @param text.size size of axis labels
#' @return ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'\dontrun{
#' p = scglmmr::GSEABarPlot(d, celltype_name = 'celltype1', fill_color = 'dodgerblue')
#' }
GSEABarPlot = function(rbind_gsea_result_dataframe, celltype_name, fill_color, text.size = 8) {
  # filter to single cell type
  dplot = rbind_gsea_result_dataframe[ rbind_gsea_result_dataframe$celltype == celltype_name, ]
  # create
  p = ggplot(dplot, aes(y = NES, x = pathway)) +
    geom_bar(stat="identity", fill=fill_color, alpha=0.8, width= 0.75, show.legend = FALSE) +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(position = "top") +
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(face = "bold",colour = "black", size = 8),
          axis.text.y = element_text(size = text.size, face = "bold", color = "black"),
          axis.text.x = element_text(size = text.size, face = "bold", color = "black")
    )
  return(p)
}






