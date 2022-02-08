# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

# these funcions are included for back compatibility with existing projects; Most are still compatible with R4 and above


###############################
# deprecated functions




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
  .Deprecated(new = 'PlotFgsea')
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


#### GSEA FUNCTIONS
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
                              nperm = 250000, positive.enrich.only = FALSE, ...) {

  .Deprecated("FgseaList")

  #require(fgsea)
  gsea = list()
  for (i in 1:length(rank.list.celltype)) {
    gsea[[i]] =
      suppressMessages(fgsea::fgsea(pathways = pathways, ...,
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







# This is a FCS repository specific function: PlotGeneDistCohort
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


  .Deprecated(c("LeadEdgeTidySampleExprs, TopGeneTidySampleExprs"))

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



