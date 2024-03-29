% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_analysis.r
\name{FgseaList}
\alias{FgseaList}
\title{FgseaList - wrapper around fast gene set enrichment analysis with the fgsea R package https://bioconductor.org/packages/release/bioc/html/fgsea.html to implement on a list of ranks indexec by cell type.}
\usage{
FgseaList(..., rank.list.celltype, pathways, maxSize = 500, minSize = 9)
}
\arguments{
\item{rank.list.celltype}{results returned by GetRankResultsRaw or GetRankResults}

\item{pathways}{modules / gene sets as a named list each a single vector of unique gene IDS}

\item{maxSize}{see fgsea package}

\item{minSize}{see fgsea package}
}
\value{
results from fgsea package indexed by celltype
}
\description{
FgseaList - wrapper around fast gene set enrichment analysis with the fgsea R package https://bioconductor.org/packages/release/bioc/html/fgsea.html to implement on a list of ranks indexec by cell type.
}
\examples{
\dontrun{
t1hvl_rank = GetRankResultsRaw(limma.fit.object.list  = dreamfit,
coefficient.number = 1,
contrast.name = "contrastName")
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
gsealist = FgseaList(rank.list.celltype = t1hvl_rank, pathways = btm,  BPPARAM = pparam)
}
# usage:
}
