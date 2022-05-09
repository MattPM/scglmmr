# scglmmr pseudobulk differential expression pipeline
# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com


#' SubjectCelltypeTable - quality control on sample level aggregation across cell types
#' @param metadata dataframe of meta data for cells-rows variables-columns i.e. ColData or
#' seurat@meta.data
#' @param sample_col quoted character e.g. "sample" the subject level sample variable - if multiple timepoints helps to code as subjectID_timepoint i.e. s1_0, s1_1
#' @param celltype_col quoted character e.g. "celltype" - the celltypes / clusters for which to create bulk libraries
#' @return a R list with table
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @examples
#'\dontrun{
#' # define counts and metadata and subset to cells above rm seurat object from workspace
#' meta = s@meta.data
#' umi = s@assays$RNA@counts
#' rm(s); gc()

#' # QC contingency of cells by subject for each celltype
#' tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype", sample_column = "sample")
#' tab$celltypes_remove; tab$`low representation celltypes`; tab$table

#' # remove cells prior to pseudobulk analysis
#' meta = meta[!meta$celltype_label_3 %in% tab$celltypes_remove, ]
#' # subset data
#' umi = umi[ ,rownames(meta)]
#' # proceed to PseudobulkList.

#' }
SubjectCelltypeTable = function(metadata, celltype_column, sample_column) {
  # make contingency table
  out = table(metadata[[celltype_column]], metadata[[sample_column]])

  # flag cell types with low or 0 representaiton
  lowcelltypes = unique(rownames(which(x = out < 10, arr.ind = TRUE, useNames = TRUE)))
  flagcelltypes = unique(rownames(which(x = out <=1, arr.ind = TRUE, useNames = TRUE)))

  # print message for user
  if(!is.null(lowcelltypes)) {
    print("see output$`low representation celltypes` with less than 10 cells for some samples:")
    message(paste0(lowcelltypes , ", " ))
  }
  if(!is.null(flagcelltypes)) {
    print('Celltypes with 1 or less cell cells for n>0 samples, remove these cells from metadata and counts prior to PseudobulkList()')
    message(paste0(flagcelltypes, ", " ))
  }

  # return list with md and draw contingency heatmap
  object = list(out , lowcelltypes, flagcelltypes)
  names(object) = c("table", "low representation celltypes", "celltypes_remove")

  # show quick heatmap
  if (!is.null(dev.list())) {
    dev.off()
    }
  pheatmap::pheatmap(t(object$table),
                     cluster_cols = FALSE, cluster_rows = FALSE,
                     fontsize_row = 6, fontsize_col = 6)
  return(object)
}

#' PseudobulkList - make a pseudobuk summed or average dataset for each samplexcelltype
#'
#' @param rawcounts the raw count UMI data for assay to sum or average for each celtype and sample.
#' @param metadata dataframe of meta data for cells-rows variables-columns i.e. ColData or seuratmeta.data
#' @param sample_col quoted character e.g. "sample" the subject level sample variable - if multiple timepoints helps to code as subjectID_timepoint i.e. s1_0, s1_1
#' @param celltype_col quoted character e.g. "celltype" - the celltypes / clusters for which to create bulk libraries
#' @param avg_or_sum whether to compute default 'sum' or 'average' library for each sample within each celltype (recommend sum)
#'
#' @return a R list of standard R matrices indexed by celltype
#' @importFrom Matrix rowSums rowMeans
#' @importFrom tibble rownames_to_column remove_rownames
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#'\dontrun{
#'pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample",
#'         celltype_col = "lineage", avg_or_sum = "sum")
#' }
PseudobulkList = function(rawcounts, metadata, sample_col, celltype_col, avg_or_sum = 'sum') {
  stopifnot(isFALSE(all(grepl(pattern = "~", unique(metadata[[celltype_col]])))))

  # make cell contingency matrix
  out = table(metadata[[celltype_col]], metadata[[sample_col]])
  flagcelltypes = unique(rownames(which(x = out < 1, arr.ind = TRUE, useNames = TRUE)))
  # warn user if celltypes without represenation in some samples have not been removed
  if (!is.null(flagcelltypes)) {
    warning('at least one sample has 0 cells for ', flagcelltypes,
            "recommend removing cells of this celltype from cell metadata and raw counts" )
  }
  # remove genes with no counts across all cells
  rawcounts = rawcounts[Matrix::rowSums(rawcounts) > 0, ]

  # list of cell barcodes indexed by sample(subject_timepoint) crossed by celltype
  metadata$celltype_sample = paste(metadata[[celltype_col]], metadata[[sample_col]], sep = "~")

  # make list of vector of cell barcodes indexed by celltype_sample
  scell = lapply(X = split(metadata, f = metadata$celltype_sample), FUN = rownames)

  # aggregate into Pseudobulk libraries as the sum or average per of cells per sample x celltype
  if (avg_or_sum == 'average') {
    csample = lapply(scell, function(x) Matrix::rowMeans(rawcounts[ ,x]))
  } else {
    csample = lapply(scell, function(x) Matrix::rowSums(rawcounts[ ,x]))
  }
  # format as a list of matrices by celltype with identical sample column names
  cmt = as.data.frame(t(do.call(cbind, csample))) %>%
    tibble::rownames_to_column("sx") %>%
    tidyr::separate("sx", into = c('celltype', 'sample'), sep = "~")
  clist = lapply(X = split(cmt, f = cmt$celltype), function(x) {
    x %>% select(-celltype) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("sample") %>%
      t()
  })
  # check column names match
  if (isTRUE(all(duplicated.default(lapply(clist, colnames)[-1])))) {
    stop(paste(' recommend removal of cell types with 0 cells for specific samples',
               ' from umi matrix and cell metadata'
    ))
  }
  return(clist)
}


#' AggregateCellMetadata - take cell level metadata and collapse down to sample-level
#' metadata for use in pseudobuk testing. This can be useful if you do not already have
#' metadata for each sample in the experiment, but these data are stored in single cell
#' metadata. For example, metadata for the donor ID, unique sample
#' (e.g. donorid_timepoint), age, sex etc. Note these are all individual sample level data.
#' Single cell intrinsic variables like nUMI cannot be collapsed down, only variables that
#' are unique for each sample which are the columns of the pseudobulk data. Only include
#' metadata that you intend to adjust models for as covariates or random effects because all
#' variables will be referenced during count normalization and feature filtering.
#'
#' @param cell.metadata dataframe of meta data for cells-rows as columns i.e. ColData or
#' Seurat@meta.data.
#' @param sample_column quoted character e.g. "sample". This should indicate the variable corresponding
#' to the rows of the pseudobulk gene expression data.
#' @param variable_columns experiment variables coresponding to sample level data. For example:
#' c('SubjectID', 'timepoint', 'sex', 'age').
#' @param pseudobulklist the object output from PseudobulkList. used to check the columns of the aggregated
#' metadata match the columns of the elements of the Pseudobulk list.
#'
#' @return Aggregated metadata
#' @importFrom dplyr group_by summarize_each
#' @export
#'
#' @examples
#'\dontrun{
#' samplemd = AggregateCellMetadata(cell.metadata = s@meta.data, sample_column = 'sample', variable_columns = c('subjectID', 'timepoint', 'age', 'sex'), pseudobulk.List = pb)
#' }
AggregateCellMetadata = function(cell.metadata, sample_column, variable_columns, pseudobulk.List){

  # define sample_column argument as data-variable to make non environment variable
  gvar = rlang::sym(sample_column)

  # collapse cells down to data.variables in model formula
  model_md = cell.metadata[ ,c(sample_column, variable_columns)] %>%
    dplyr::group_by({{gvar}}) %>%
    dplyr::summarise_each(list(~unique(.)))

  # convert tibble to dataframe to use in dream fit
  model_md = base::as.data.frame(model_md)
  rownames(model_md) = model_md[[sample_column]]

  # quality Control against the pseudobulk list per cell type
  # rows match the pseudobulk data columns
  expr.colnames = colnames(pseudobulklist[[1]])
  if (isFALSE(all.equal(target = expr.colnames, current = rownames(model_md)))) {
    stop('rows of design matrix do not match column names of pseudobulklist')
  }
  # show and return design matrix
  print(model_md); print(str(model_md))
  return(met)
}


#' Normalize - normalize summed gene counts across samples. This is a wrapper around edgeR functions `calcNormFactors` and `filterByExpr`.
#' @param pseudobulk.list object created with PseudobulkList only use this function if argument to PseudobulkList avg_or_sum was 'sum' computes normalization for pseudobulk libraries
#' @param normalization.method see edgeR function calcNormFactors, this is the argument to `method`
#' @param design see edgeR `filterByExpr` function. This is the design matrix created using model.matrix and e.g. the metadata created  AggregateCellMetadata.
#' @param minimum.gene.count see edgeR function `filterbyExpr` thie is the argument to `min.count`
#'
#' @return a list of dgeList objects indexed by celltype
#' @importFrom edgeR calcNormFactors filterByExpr
#' @export
#'
#' @examples
#' #' #'\dontrun{
# 'dge = scglmmr::Normalizek(pseudobulk.list = pb, sample.metadata = metadata, minimum.gene.count = 5, normalization.method = 'TMM')
#' }
Normalize = function(pseudobulk.list,
                     design,
                     group = NULL,
                     normalization.method = "RLE",
                     minimum.gene.count = 1) {
  # convert object to DGEList
  dflist = lapply(pseudobulklist, edgeR::DGEList, samples = sample.metadata)

  # calculate normalization factors
  dflist = lapply(dflist, function(x)
    edgeR::calcNormFactors(x, method = normalization.method))

  # define expressed features for each cell type
  genes.use1 = lapply(dflist, function(x) {
    edgeR::filterByExpr(x, min.count = minimum.gene.count, design = design)
  })

  # filter features that are not expressed.
  for (i in 1:length(pseudobulklist)) {
    print(names(pseudobulklist)[i])
    print('retained genes after filtering:  ')
    print(nrow(dflist[[i]]$counts))
    # filter the dge lists and recalculate library size
    dflist[[i]] = dflist[[i]][genes.use1[[i]], keep.lib.sizes = FALSE]
  }
  return(dflist)
}




###############################
# Deprecated functions

#' BulkDesignMatrix - designmatrix for the main experiment factor - used in normalization
#'
#' @param metadata dataframe of meta data for cells-rows variables-columns i.e. ColData or seuratmeta.data
#' @param sample_column quoted character e.g. "sample" the subject level sample variable - if multiple timepoints should be subjectID_timepoint i.e. s1_0, s1_1
#' @param variable_column the main experiment variable of interest e.g. timepoint or if implementing custom contrasts for difference in foldchange between groups, a combined factor of group_time, e.g. pooroutcome_t0, pooroutcome_t1, goodoutcome_t0, goodoutcome_t1 additional covariates can be specified in dreamMixedModel
#' @param pseudobulklist the list created with PseudobulkList
#'
#' @return a design matrix
#' @import Matrix
#' @export
#'
#' @examples
#'\dontrun{
#' designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",
#' variable_column = "cohort_timepoint", pseudobulklist = pb)
#' }
BulkDesignMatrix = function(metadata, sample_column, variable_column, pseudobulklist){
  .Deprecated(new = 'AggregateCellMetadata')
  # create design matrix from aggregated contingency table of cells
  met = as.data.frame.matrix(table(metadata[[sample_column]], metadata[[variable_column]]))
  # convert contingency count to 1 for design matrix
  met[met>0] = 1

  # add labels for variancepartition contrast fit.
  colnames(met) = paste(variable_column,colnames(met), sep = "")

  # QC design matrix
  # rows match the pseudobulk data columns
  if (isFALSE(all.equal(target = colnames(pseudobulklist[[1]]), current = rownames(met)))) {
    stop('rows do not match column names of pseudobulklist')
  }
  # show and return design matrix
  met
  return(met)
}

#' NormalizePseudobulk
#'
#' @param pseudobulklist object created with PseudobulkList only use this function if argument to PseudobulkList avg_or_sum was 'sum' computes normalization for pseudobulk libraries
#' @param normalization.method see edgeR function calcNormFactors, this is the argument to `method`
#' @param design_matrix the sample level metadata created by AggregateCellMetadata which will be used to make the designamatrix.
#' @param minimum.gene.count see edgeR function `filterbyExpr` thie is the argument to `min.count`
#'
#' @return a list of dgeList indexed by celltype
#' @importFrom edgeR calcNormFactors filterByExpr
#' @export
#'
#' @examples
#' #' #'\dontrun{
# 'dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)
#' }
#' # normalize the pseudobulk data made with MakePseudobulkList
NormalizePseudobulk = function(pseudobulklist,
                               normalization.method = "RLE",
                               design_matrix,
                               minimum.gene.count = 1) {
  .Deprecated('Normalize')
  dflist = lapply(pseudobulklist, edgeR::DGEList)
  dflist = lapply(dflist, function(x)
    edgeR::calcNormFactors(x, method = normalization.method))
  genes.use1 = lapply(dflist, function(x) {
    edgeR::filterByExpr(x, min.count = minimum.gene.count, design = design_matrix)
  })
  for (i in 1:length(pseudobulklist)) {
    print(names(pseudobulklist)[i])
    print('retained genes after filtering:  ')
    print(nrow(dflist[[i]]$counts))
    dflist[[i]] = dflist[[i]][genes.use1[[i]], keep.lib.sizes = FALSE]
  }
  return(dflist)
}






