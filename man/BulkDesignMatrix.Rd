% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_helpers.r
\name{BulkDesignMatrix}
\alias{BulkDesignMatrix}
\title{BulkDesignMatrix - designmatrix for the main experiment factor - used in normalization}
\usage{
BulkDesignMatrix(metadata, sample_column, variable_column, pseudobulklist)
}
\arguments{
\item{metadata}{dataframe of meta data for cells-rows variables-columns i.e. ColData or seuratmeta.data}

\item{sample_column}{quoted character e.g. "sample" the subject level sample variable - if multiple timepoints should be subjectID_timepoint i.e. s1_0, s1_1}

\item{variable_column}{the main experiment variable of interest e.g. timepoint or if implementing custom contrasts for difference in foldchange between groups, a combined factor of group_time, e.g. pooroutcome_t0, pooroutcome_t1, goodoutcome_t0, goodoutcome_t1 additional covariates can be specified in dreamMixedModel}

\item{pseudobulklist}{the list created with PseudobulkList}
}
\value{
a design matrix
}
\description{
BulkDesignMatrix - designmatrix for the main experiment factor - used in normalization
}
\examples{
\dontrun{
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",
variable_column = "cohort_timepoint", pseudobulklist = pb)
}
}
