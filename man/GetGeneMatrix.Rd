% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_result_interaction.r
\name{GetGeneMatrix}
\alias{GetGeneMatrix}
\title{GetGeneMatrix get a gene matric for plotting genes by celltypes statistic from pseudobulk model results e.g. using heatmap pheatmap or complexheatmap}
\usage{
GetGeneMatrix(
  result.list,
  gene_subset = NULL,
  stat_for_matrix = "logFC",
  pvalfilter,
  logfcfilter
)
}
\arguments{
\item{result.list}{the object returned by GetContrastResults or GetContrastResultsRaw}

\item{gene_subset}{a preselected subset of genes as a vector, for example returned by GetLeadingEdgeFull. Defaults to the union of all fitted genes across cell types.}

\item{stat_for_matrix}{defaults to logFC, the effect size. can be any of the columns returned by limma::topTable}

\item{pvalfilter}{filter genes to retain in the matrix by the raw p values}

\item{logfcfilter}{filter genes to retain in the matrix by the logFC}
}
\value{
a matrix
}
\description{
GetGeneMatrix get a gene matric for plotting genes by celltypes statistic from pseudobulk model results e.g. using heatmap pheatmap or complexheatmap
}
\examples{
\dontrun{
le = scglmmr::GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.1, NES.filter = -Inf)
genesub = do.call(rbind, le) \%$\% gene \%>\% unique
mtx2 = scglmmr::GetGeneMatrix(result.list = res,
                              stat_for_matrix = "logFC",
                              gene_subset = genesub,
                              pvalfilter = -Inf,
                              logfcfilter = 0.1)

pheatmap::pheatmap(mtx2,
                   breaks =seq(from = 0, to = 2,length.out = 99),
                   filename = paste0(figpath,"LEgenes_heatmap.pdf"))
}
}
