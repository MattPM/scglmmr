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

