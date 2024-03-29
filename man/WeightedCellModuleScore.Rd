% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_cell_gene_module_scores.r
\name{WeightedCellModuleScore}
\alias{WeightedCellModuleScore}
\title{WeightedCellModuleScore - calculate the average or weighted aveage +/- scaling return dataframe of cells by average module score}
\usage{
WeightedCellModuleScore(
  gene_matrix = NULL,
  module_list,
  threshold = 0,
  cellwise_scaling = FALSE,
  return_weighted = FALSE
)
}
\arguments{
\item{gene_matrix}{- normalized genes (rows) by cells (columns)}

\item{module_list}{- names list of gene modules - each element a vector of gene names.}

\item{threshold}{- at least this fraction of genes in the signature must be < 0 across all cells or else a score is not calculated}

\item{cellwise_scaling}{- T/F scale across cells ?}

\item{return_weighted}{- T/F weight the averae by multiplying by gene representation?}
}
\value{
dataframe of cells barcodes (rownames) by gene module scores (columns)
}
\description{
WeightedCellModuleScore - calculate the average or weighted aveage +/- scaling return dataframe of cells by average module score
}
\examples{
\dontrun{
}
}
