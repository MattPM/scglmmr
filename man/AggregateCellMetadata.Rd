% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_helpers.r
\name{AggregateCellMetadata}
\alias{AggregateCellMetadata}
\title{AggregateCellMetadata - take cell level metadata and collapse down to sample-level
metadata for use in pseudobuk testing. This can be useful if you do not already have
metadata for each sample in the experiment, but these data are stored in single cell
metadata. For example, metadata for the donor ID, unique sample
(e.g. donorid_timepoint), age, sex etc. Note these are all individual sample level data.
Single cell intrinsic variables like nUMI cannot be collapsed down, only variables that
are unique for each sample which are the columns of the pseudobulk data. Only include
metadata that you intend to adjust models for as covariates or random effects because all
variables will be referenced during count normalization and feature filtering.}
\usage{
AggregateCellMetadata(
  cell.metadata,
  sample_column,
  variable_columns,
  pseudobulk.List
)
}
\arguments{
\item{cell.metadata}{dataframe of meta data for cells-rows as columns i.e. ColData or
Seurat@meta.data.}

\item{sample_column}{quoted character e.g. "sample". This should indicate the variable corresponding
to the rows of the pseudobulk gene expression data.}

\item{variable_columns}{experiment variables coresponding to sample level data. For example:
c('SubjectID', 'timepoint', 'sex', 'age').}

\item{pseudobulklist}{the object output from PseudobulkList. used to check the columns of the aggregated
metadata match the columns of the elements of the Pseudobulk list.}
}
\value{
Aggregated metadata
}
\description{
AggregateCellMetadata - take cell level metadata and collapse down to sample-level
metadata for use in pseudobuk testing. This can be useful if you do not already have
metadata for each sample in the experiment, but these data are stored in single cell
metadata. For example, metadata for the donor ID, unique sample
(e.g. donorid_timepoint), age, sex etc. Note these are all individual sample level data.
Single cell intrinsic variables like nUMI cannot be collapsed down, only variables that
are unique for each sample which are the columns of the pseudobulk data. Only include
metadata that you intend to adjust models for as covariates or random effects because all
variables will be referenced during count normalization and feature filtering.
}
\examples{
\dontrun{
samplemd = AggregateCellMetadata(cell.metadata = s@meta.data, sample_column = 'sample', variable_columns = c('subjectID', 'timepoint', 'age', 'sex'), pseudobulk.List = pb)
}
}
