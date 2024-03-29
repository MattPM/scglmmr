% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_deprecated_functions.r
\name{GSEABubblePlot}
\alias{GSEABubblePlot}
\title{GSEABubblePlot plot gsea results for all cell types}
\usage{
GSEABubblePlot(
  rbind_gsea_result_dataframe,
  save_path,
  include_negative = TRUE,
  save_name,
  width = 8.5,
  height = 7.2
)
}
\arguments{
\item{rbind_gsea_result_dataframe}{dataframe returned by RbindGseaResultList}

\item{save_path}{file path to save results}

\item{include_negative}{TRUE/FALSE whether to include negative enrichment in the plot.}

\item{save_name}{name of file saved to save_path}

\item{width}{ggpsave param}

\item{height}{ggsave param}
}
\value{
nothing
}
\description{
GSEABubblePlot plot gsea results for all cell types
}
\examples{
\dontrun{
scglmmr::GSEABubblePlot(d, save_path = figpath, save_name = "plot.pdf")
}
}
