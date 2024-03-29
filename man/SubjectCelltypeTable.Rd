% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk_helpers.r
\name{SubjectCelltypeTable}
\alias{SubjectCelltypeTable}
\title{SubjectCelltypeTable - quality control on sample level aggregation across cell types}
\usage{
SubjectCelltypeTable(metadata, celltype_column, sample_column)
}
\arguments{
\item{metadata}{dataframe of meta data for cells-rows variables-columns i.e. ColData or
seurat@meta.data}

\item{sample_col}{quoted character e.g. "sample" the subject level sample variable - if multiple timepoints helps to code as subjectID_timepoint i.e. s1_0, s1_1}

\item{celltype_col}{quoted character e.g. "celltype" - the celltypes / clusters for which to create bulk libraries}
}
\value{
a R list with table
}
\description{
SubjectCelltypeTable - quality control on sample level aggregation across cell types
}
\examples{
\dontrun{
# define counts and metadata and subset to cells above rm seurat object from workspace
meta = s@meta.data
umi = s@assays$RNA@counts
rm(s); gc()
# QC contingency of cells by subject for each celltype
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype", sample_column = "sample")
tab$celltypes_remove; tab$`low representation celltypes`; tab$table
# remove cells prior to pseudobulk analysis
meta = meta[!meta$celltype_label_3 \%in\% tab$celltypes_remove, ]
# subset data
umi = umi[ ,rownames(meta)]
# proceed to PseudobulkList.
}
}
