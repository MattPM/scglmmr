% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_deprecated_functions.r
\name{GetContrastResultsRaw}
\alias{GetContrastResultsRaw}
\title{GetContrastResultsRaw - calculate p values and return contrast results from modelfit with dreamMixedModel}
\usage{
GetContrastResultsRaw(limma.fit.object.list, coefficient.number, contrast.name)
}
\arguments{
\item{limma.fit.object.list}{the results returned by dreamMixedModel, to get coefficiennt from RunVoomLimma use GetContrastResults}

\item{coefficient.number}{corresponds to the contrast, the nmber is in order of the contrast matrix}

\item{contrast.name}{this does not have to match the exact contrast name used in the contrast matrix, it is here to force / remind the user to choose the correct contrast.}
}
\value{
a list of dataframes with contrast results indexed by cell type
}
\description{
GetContrastResultsRaw - calculate p values and return contrast results from modelfit with dreamMixedModel
}
\examples{
\dontrun{
fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit,
                                        coefficient.number = 1,
                                        contrast.name = "foldchangedifference")

}
}
