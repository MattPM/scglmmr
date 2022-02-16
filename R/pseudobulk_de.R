# scglmmr pseudobulk differential expression pipeline
# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

#' SubjectCelltypeTable - make a QC contingency table of cells by samples; this function sets up the workflow and is run prior to running PseudobulkList
#'
#' @param metadata metadata cells-rows variables-columns i.e. ColData or seurat meta.data
#' @param celltype_column quoted character e.g. "celltype" - the celltypes / clusters for which to create bulk libraries
#' @param sample_column  quoted character e.g. "sample" the subject level sample variable - if multiple timepoints should be subjectID_timepoint i.e. s1_0, s1_1
#'
#' @return a list with index 1 = contingency table of cells by sample, 2 = low representation celltypes < 10 cels, 3 celltypes to remove with 0 cells for some samples.
#' @importFrom stats heatmap
#' @export
#'
#' @examples
#'\dontrun{
#' datapath = "mypath/"
#'
#' # load seurat or sce object etc.
#' s = readRDS("path/seuratobject.rds")
#'
#' # define counts and metadata and subset to cells above rm seurat object from workspace
#' # meta = s@meta.data
#' # umi = s@assays$RNA@counts
#' # rm(s); gc()
#'
#' # QC contingency of cells by subject for each celltype
#' tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "lineage", sample_column = "sample")
#' tab$celltypes_remove; tab$`low representation celltypes`; tab$table
#'
#' # remove cells prior to pseudobulk analysis
#' meta = meta[!meta$celltype_label_3 %in% tab$celltypes_remove, ]
#'
#' # subset data
#' umi = umi[ ,rownames(meta)]
#' # proceed to PseudobulkList
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
  stats::heatmap(t(object$table), Rowv = NA, Colv = NA, keep.dendro = FALSE)
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
#' @import Matrix
#' @importFrom tibble rownames_to_column remove_rownames
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#'\dontrun{
#'pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample",
#'         celltype_col = "lineage", avg_or_sum = "sum")
#' }
#' # for each celltype aggregate raw counts for ech sample (subject timepoint combination)
# store in a list of matrices
PseudobulkList = function(rawcounts, metadata, sample_col, celltype_col, avg_or_sum = 'sum') {
  stopifnot(isFALSE(all(grepl(pattern = "~", unique(metadata[[celltype_col]])))))

  # make cell contingency matrix
  out = table(metadata[[celltype_col]], metadata[[sample_col]])
  flagcelltypes = unique(rownames(which(x = out < 1, arr.ind = TRUE, useNames = TRUE)))
  # warn user if celltypes without represenation in some samples have not been removed
  if (!is.null(flagcelltypes)) {
    warning('some samples have 0 cells for ', flagcelltypes,
            "cannot use RunVoomLimma or dreamMixedModel without custom design matrix for each celltype",
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
    stop(paste('samples not identical for all celltypes.',
               ' remove celltypes with 0 cells in some samples from counts and metadata'
               ))
  }
  return(clist)
}


#' BulkDesignMatrix - designmatrix for the main experiment factor
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

  # create design matrix from aggregated contingency table of cells
  met = as.data.frame.matrix(table(metadata[[sample_column]], metadata[[variable_column]]))
  # convert contingency count to 1 for design matrix
  met[met>0] = 1

  # add labels for variancepartition contrast fit.
  colnames(met) = paste(variable_column,colnames(met), sep = "")

  # QC design matrix
  # rows match the pseudobulk data columns
  if (isFALSE(all.equal(target = colnames(pseudobulklist[[1]]), current = rownames(met)))) {
    stop('rows of design matrix do not match column names of pseudobulklist')
  }

  # check if the design matrix is full rank
  if(isFALSE(Matrix::rankMatrix(met) == ncol(met))){
    warning('design matrix is not full rank')
  }

  # show and return design matrix
  met
  return(met)
}


#' NormalizePseudobulk
#'
#' @param pseudobulklist object created with PseudobulkList only use this function if argument to PseudobulkList avg_or_sum was 'sum' computes normalization for pseudobulk libraries
#' @param normalization.method see edgeR function calcNormFactors, this is the argument to `method`
#' @param design_matrix the design matrix created with `BulkDesignMatrix`
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
NormalizePseudobulk = function(pseudobulklist, normalization.method = "RLE",
                               design_matrix, minimum.gene.count = 1) {
  require(edgeR)
  dflist = lapply(pseudobulklist, edgeR::DGEList)
  dflist = lapply(dflist, function(x) edgeR::calcNormFactors(x, method = normalization.method))
  genes.use1 = lapply(dflist, function(x) {
    edgeR::filterByExpr(x, min.count = minimum.gene.count, design = design_matrix)
  })
  for (i in 1:length(pseudobulklist)) {
  print(names(pseudobulklist)[i])
  print('retained genes after filtering:  ')
  dflist[[i]] = dflist[[i]][genes.use1[[i]], keep.lib.sizes=FALSE]
  print(nrow(dflist[[i]]$counts))
  }
  return(dflist)
}


#' RunVoomLimma
#'
#' @param dgelists a ist of DGElist created in NormalizePseudobulk
#' @param design_matrix design matrix created with BulkDesignMatrix
#' @param do_contrast_fit whether to fit custom a priori contrasts
#' @param my_contrast_matrix custom a priori contrasts created with make.contrasts - see limma or edgeR manual
#'
#' @return a list of linear model fits for each celltype
#' @importFrom limma lmFit voom contrasts.fit eBayes
#' @export
#'
#' @examples
#'\dontrun{
#'c_mat = makeContrasts(
#'foldchange_difference = (group_timepoint1_1 - group_timepoint1_0) - (group_timepoint0_1 - group_timepoint0_0),
#'time1_foldchange = (group_timepoint1_1 + group_timepoint0_1) / 2  - (group_timepoint1_0 + group_timepoint0_0) / 2,
#'baseline_groups = (group_timepoint1_0 - group_timepoint0_0),
#'levels = colnames(designmat)
#')
#'# fit simple linear model for the baseline group level contrast
#'bl = scglmmr::RunVoomLimma(dgelists = dge,
#'design_matrix = designmat,
#'do_contrast_fit = T,
#'my_contrast_matrix = c_mat[ ,3])
#' }
#' # run limma using voom observational weights for non mixed effects models using emperical bayes
RunVoomLimma = function(dgelists, design_matrix, do_contrast_fit, my_contrast_matrix){
  print("to implement random intercept for repeated measures (e.g. time) from same donor use dreamMixedModel")

  # get voom observational weights
  v = lapply(dgelists, function(x){
    limma::voom(counts = x, design = design_matrix,
                normalize.method = "none",save.plot = T, plot = F)
    })
  # fit model and apply custom contrasts if specified
  fit = lapply(v, function(x){
    limma::lmFit(x, design = design_matrix)
    })
    if (do_contrast_fit == TRUE) {
      c_fit = lapply(fit, function(x){
        limma::contrasts.fit(x, contrasts = my_contrast_matrix)
        })
      eb_c_fit = lapply(c_fit, limma::eBayes)
    } else {
      eb_c_fit = lapply(fit, limma::eBayes)
    }
  # name result list by celltype
  names(eb_c_fit) = names(dgelists)
  return(eb_c_fit)
}


#' dreamMixedModel - run dream mixed model
#' @param dge_lists list of dgelists created with NormalizePseudobulk
#' @param apriori_contrasts one of TRUE or FALSE, whether to fit a priori contrasts

#' @param contrast_matrix contrast matrix created with make.contrasts
#' @param design_matrix design matrix created with BulkDesignMatrix
#' @param plotsavepath a path to save created plot of contrasts
#' @param ncores number of cores for doParallel
#' @param cell_metadata metadata cells-rows variables-columns i.e. ColData or seuratmeta.data
#' @param fixed_effects a vector of covariates that are columns in metadata e.g. a vector of  'age' 'gender'
#' @param sample_column quoted character e.g. "sample" the subject level sample variable should have multiple timepoints subjectID_timepoint i.e. s1_0, s1_1
#' @param lme4_formula symbolic model formula the default is '~ 0 + cohort_timepoint + (1|sampleid)'
#'
#' @return list of model fits indexed by celltype
#' @importFrom variancePartition voomWithDreamWeights plotContrasts dream
#' @importFrom limma voom
#' @importFrom parallel makeCluster
#' @importFrom ggplot2 theme ggsave
#' @importFrom BiocParallel register SnowParam
#' @importFrom rlang sym
#' @importFrom dplyr group_by summarize_each
#' @export
#'
#' @examples
#'\dontrun{
#'c_mat = makeContrasts(
#'foldchange_difference = (group_timepoint1_1 - group_timepoint1_0) - (group_timepoint0_1 - group_timepoint0_0),
#'time1_foldchange = (group_timepoint1_1 + group_timepoint0_1) / 2  - (group_timepoint1_0 + group_timepoint0_0) / 2,
#'baseline_groups = (group_timepoint1_0 - group_timepoint0_0),
#'levels = colnames(designmat)
#')
#'# fit mixed model for the multi timepoint contrasts
#'fit = scglmmr::dreamMixedModel(dge_lists = dge,
#'apriori_contrasts = TRUE,
#'sample_column = 'sample',
#'cell_metadata = meta,
#'contrast_matrix = c_mat,
#'design_matrix = designmat,
#'lme4_formula =  '~ 0 + age + gender + cohort_timepoint + (1|sampleid)',
#'fixed_effects = c('age', 'gender', 'cohort_timepoint'),
#'plotsavepath = figpath,
#'ncores = 4)
#' }
#' # run dream mixed model
# dream method Hoffman et.al. 2020  https://doi.org/10.1093/bioinformatics/btaa687
dreamMixedModel = function(dge_lists, apriori_contrasts = FALSE, sample_column, contrast_matrix = NULL,
                           design_matrix, fixed_effects, cell_metadata,
                           lme4_formula = '~ 0 + cohort_timepoint + (1|sampleid)',
                           ncores = 4, ...) {

  BiocParallel::register(BiocParallel::SnowParam(workers = ncores))
  pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
  # cl = parallel::makeCluster(ncores)
  # doParallel::registerDoParallel(cl = ncores)

  if(isFALSE('sampleid' %in% colnames(cell_metadata))){
    stop("metadata must have a column 'sampleid' to fit a varying intercept model with subjectIDs")
  }

  # combine vector of covariates in model
  vars_all = c(sample_column, fixed_effects, 'sampleid')

  # define sample_column argument as data-variable to make non environment variable
  gvar = rlang::sym(sample_column)
  # collapse cells down to data.variables in model formula {{}} for non standard evaluation of data.variable
  model_md = cell_metadata[ ,vars_all] %>%
    dplyr::group_by({{gvar}}) %>%
    dplyr::summarise_each(list(~unique(.)))

  # convert tibble to dataframe to use in dream fit
  model_md = base::as.data.frame(model_md)
  rownames(model_md) = model_md[[sample_column]]
  print("dream data argument for model "); print(model_md)
  # check row index of model metadata matches dgelist columns and show symbolic formula
  stopifnot(isTRUE(all.equal(target = colnames(dge_lists[[1]]), current = rownames(model_md))))

  print('model specified (change with argument to lme4_formula) '); print(lme4_formula)
  v1 = lapply(dge_lists, function(x){
    design = NULL
    variancePartition::voomWithDreamWeights(counts = x, data = model_md, formula = lme4_formula,
                                            normalize.method = "none", save.plot = TRUE, plot = TRUE)
    })
  if(isTRUE(apriori_contrasts)) {
    # fit mixed model
    fit1 = lapply(v1, function(x){
      variancePartition::dream(...,
                               exprObj = x,
                               formula = lme4_formula,
                               data = model_md,
                               L = as.matrix(contrast_matrix),
                               BPPARAM = pparam)
    })
  } else {
    # fit model without custom contrasts
    fit1 = lapply(v1, function(x){
      variancePartition::dream(exprObj = x, formula = lme4_formula, data = model_md)
    })
  }
  return(fit1)
}

#' calc_avg_module_zscore calculate average module z score of list of modules on a PseudobulkList
#'
#' @param module.list list of modules
#' @param average.data.frame - this is created in AverageSampleModuleZscore
#'
#' @return see AverageSampleModuleZscore
#' @export
#'
#' @examples
#'\dontrun{
#' results = calc_avg_module_zscore(module.list = btm, average.data.frame = av_df)
#' }
#' # Average Module sample Z score
# the method below is equivalent to Yuri's function used in baseline paper Kotliarov et. al. Nat Med 2020
# it is adopted below to run on 'pseudobulk lists' (average "averagemetacell.list" or pseudobulk list created by PseudobulkList)
# The calc_avg_module_zscore function gets called in the AverageSampleModuleZscore function below.
# calculate signature score for each cell type, BTM, Subject
# function input = named list of modules, dataframe with subject as rows genes as columns
calc_avg_module_zscore = function(module.list, average.data.frame) {
  res = data.frame()
  for (u in 1:length(module.list)) {
    # subset data by genes in module
    av = average.data.frame %>% base::as.matrix()
    mod.genes = module.list[[u]] %>% as.vector
    mod.genes = rownames(av) %in% mod.genes
    av = av[mod.genes, ]

    # scale genes for each subject, get average of z score
    x = av %>% t() %>% scale() %>% t()
    x = colMeans(x, na.rm=T) %>% t() %>% as.data.frame()
    res = rbind(res, x)
  }
  rownames(res) = names(module.list)
  return(res)
}


#' AverageSampleModuleZscore apply the function  calc_avg_module_zscore to a pseudobulklist.
#'
#' @param average.metacell.list poorly named argument - the object created by PseudobulkList either an average or summed pseudobulk data
#' @param module.list list of modules as named list each element is a vector of gene names
#' @param use.module.subset TRUE or FALSE - calc a different set of modules for each celltype use with modules.subset.by.celltype
#' @param modules.subset.by.celltype modules.subset.by.celltype is a list of modules to test with length = celltypes.vector and n modules = unique to subset
#'
#' @return returns a dataframe of module scores for each celltype
#' @export
#'
#' @examples
#'\dontrun{
#'AverageSampleModuleZscore = function(
#'  average.metacell.list = av,
#'  module.list = btm ,
#'  use.module.subset = FALSE,
#'  )
#'
#' }
AverageSampleModuleZscore = function(average.metacell.list,
                                     module.list,
                                     use.module.subset = TRUE,
                                     modules.subset.by.celltype = modules.test) {
  mod.scores.celltype = list()
  if(use.module.subset == TRUE){
    for (i in 1:length(average.metacell.list)) {
      mod.scores.celltype[[i]] =
        calc_avg_module_zscore(module.list = module.list[modules.subset.by.celltype[[i]]],
                               average.data.frame = average.metacell.list[[i]])
    }
  } else {
    for (i in 1:length(average.metacell.list)) {
      mod.scores.celltype[[i]] =
        calc_avg_module_zscore(module.list = module.list[modules.subset.by.celltype[[i]]],
                               average.data.frame = average.metacell.list[[i]])
    }
  }
  names(mod.scores.celltype) = celltypes
  return(mod.scores.celltype)
}


################# AJM function

# RunVoomLimma_covar - Note:
# removed celltypes vector argument to match updated function above; argument was not needed.
# same as above [RunLimmaVoom] with option to adjust for covariate within the mRNA (AJM)


#' RunVoomLimma_covar RunLimmaVoom with an option to adjust for a covariate within mRNA.
#'
#' @param dgelists DGElists created with PseudobulkList
#' @param design_matrix design matrix created with BulkDesignMatrix
#' @param co_variable_genes co variable genes
#' @param grouptable group table
#' @param grouptable argument to gsva
#' @param do_contrast_fit do contrasts
#' @param my_contrast_matrix contrast matrix
#' @param my_model_metadata model metadata
#' @param celltypes.vector vector of celltypes assign as names of dgeList
#'
#' @return
#' @importFrom GSVA gsva
#' @importFrom limma lmFit contrasts.fit eBayes voom
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#'\dontrun{
#' RunVoomLimma_covar(dgelists = dge_lists, design_matrix = design, co_variable_genes = NULL, grouptable = grouptable ,do_contrast_fit = TRUE,
#'  my_contrast_matrix = c_mat, my_model_metadata = md, celltypes.vector = NULL, parallel.sz = 4)
#'
#' }
RunVoomLimma_covar = function(dgelists, design_matrix, co_variable_genes = NULL, grouptable,
                        do_contrast_fit, my_contrast_matrix, my_model_metadata, celltypes.vector = NULL, parallel.sz = 4){

  v = cor = fit = c_fit = eb_c_fit = dt = gt = dm2 = list()
  for (i in 1:length(dgelists)) {
      # voom call 1 get logcpm observational weights to feed into dup cor
      v[[i]] <- limma::voom(counts = dgelists[[i]], normalize.method = "none", design = design_matrix, save.plot = T, plot = F)

      # covariable genes GSVA
      cor[[i]] = as.numeric(GSVA::gsva(v[[i]]$E, list(covar = co_variable_genes), parallel.sz = parallel.sz))

      gt[[i]] <- grouptable

      gt[[i]]$gsva <- cor[[i]]

      # make designmatrix2
      dm2[[i]] = stats::model.matrix(~0+gsva+group, data=gt[[i]])

      # Fit model and contrasts
      fit[[i]] = limma::lmFit(v[[i]], design = dm2[[i]])

      if (do_contrast_fit == TRUE) {
      c_fit[[i]] = limma::contrasts.fit(fit = fit[[i]], contrasts = my_contrast_matrix)
      eb_c_fit[[i]] = limma::eBayes(c_fit[[i]])
    } else {
      eb_c_fit[[i]] = limma::eBayes(fit[[i]])
    }
  }
  # if (is.null(celltypes.vector)) {
  #   return(eb_c_fit)
  # } else {
  #   names(eb_c_fit) = celltypes.vector
  #
  # }
return(eb_c_fit)
}




############## functions with seurat versioning dependency
# AJM function for metacells
# MakeMetaCellperSampleList = function(seurat.object, cluster.Res, sample_column,
#                                      num.pcs=20, assays, CITEprots, vector.of.samples,
#                                      numthreads = 1) {
#   sl = list()
#   distlist = list()
#   MetaCelllist = list()
#   samples = vector.of.samples
#   Idents(seurat.object) <- sample_column
#   for (i in 1:length(samples)) {
#     sl[[i]] = seurat.object %>%
#       subset(idents = samples[i]) %>%
#       FindVariableFeatures(nfeatures=2000) %>%
#       ScaleData(assay="RNA") %>%
#       RunPCA(assay="RNA", slot = "scale.data")
#     distlist[[i]] <- parDist(t(rbind(GetAssayData(sl[[i]][["CITE"]])[CITEprots,],t(Embeddings(sl[[i]], reduction="pca")[,num.pcs]))), threads = numthreads)
#
#     sl[[i]][["FiltCATres_snn"]] <- FindNeighbors(distlist[[i]])$snn
#     sl[[i]] <- FindClusters(sl[[i]], resolution = cluster.Res, graph.name = "FiltCATres_snn", algorithm = 1)
#   }
#   for (j in 1:length(samples)) {
#     MetaCelllist[[j]] = AverageExpression3(sl[[j]], assays = assays, features = NULL, return.seurat = TRUE,
#                                            slot = "counts", use.counts=TRUE, verbose = TRUE)
#   }
#   names(MetaCelllist) = samples
#   return(MetaCelllist)
# }


##################################
# RETIRED FUNCTIONS             #
#################################


####################
# dreamMixedModel(
# dge_lists = dflist,
# apriori_contrasts = TRUE,
# contrast_matrix = contrast_matrix,
# design_matrix =  met,
# plotsavepath = here(""),
# ncores = 4,
# cell_metadata = metadata,
# sample_column = 'sample',
# fixed_effects = c('cohort_timepoint'),
# lme4_formula = '~ 0 + cohort_timepoint + (1|sampleid)'
# )
####################


# get contrast resuls
# GetContrastResults = function(limma.fit.object.list, coefficient.number, contrast.name, celltypes.vector = celltypes){
#   require(limma)
#   test = list()
#   for (i in 1:length(limma.fit.object.list)) {
#     test[[i]] =
#       topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf) %>%
#       rownames_to_column("gene") %>%
#       mutate(contrast = rep(contrast.name)) %>%
#       mutate(pvalue0.01 = if_else(P.Value < 0.01, true = "1", false = "0")) %>%
#       mutate(celltype = celltypes.vector[i])
#   }
#   names(test) = celltypes.vector
#   return(test)
# }


# ## make a table of top DE genes acros all celltypes
#   MakeGeneTable = function(result.list, celltypes.vector, pvalthresh, logfcgreater) {
#     merged = data.frame()
#     for (i in 1:length(result.list)) {
#     	sub =
#     	result.list[[i]] %>%
#     	mutate(celltype = celltypes.vector[i]) %>%
#     	filter(P.Value  < pvalthresh) %>%
#     	filter(logFC > logfcgreater)
#       merged = rbind(merged, sub)
#     }
#     return(merged)
#   }
#
#
#
# # add number of genes parameter.
# GetTopGenesFDR = function(result.list, celltypes.vector, pvalthresh, logfcgreater, ngenes) {
#   merged = data.frame()
#   for (i in 1:length(result.list)) {
#     sub =
#       result.list[[i]] %>%
#       mutate(celltype = celltypes.vector[i]) %>%
#       filter(adj.P.Val < pvalthresh) %>%
#       filter(logFC > logfcgreater)
#     sub = sub[1:ngenes, ]
#     merged = rbind(merged, sub)
#   }
#
#   df.gene =
#     merged %>% group_by(gene) %>%
#     summarise(min.p = min(adj.P.Val, na.rm=T)) %>%
#     ungroup() %>%
#     filter(min.p < pvalthresh)
#   df.res =
#     merged %>%
#     filter(gene %in% df.gene$gene)
#   pl =
#     df.res %>%
#     select(gene, logFC, celltype) %>%
#     spread(key = celltype, value = logFC)
#   plot =
#     pl %>%
#     column_to_rownames("gene") %>%
#     as.matrix()
#   plot[is.na(plot)] = 0
#   return(plot)
# }



#
# # get rank results based on t statistic for GSEA.
# GetRankResults = function(limma.fit.object.list, coefficient.number, contrast.name){
#   require(limma)
#   ranks = list()
#   for (i in 1:length(limma.fit.object.list)) {
#   test = as.data.frame(topTable(limma.fit.object.list[[i]], coef = coefficient.number, number = Inf))
#   ranks[[i]] =
#     test %>% rownames_to_column("gene") %>%
#     arrange(desc(t)) %>%
#     select(gene, t) %>%
#     column_to_rownames("gene") %>%
#     t() %>%
#     unlist(use.names = T)
#   ranks[[i]] = ranks[[i]][1, ]
#   }
#   return(ranks)
# }
#
#
# # run GSEA wrapper around fgsea function
# RunFgseaOnRankList = function(rank.list.celltype, pathways, celltypes.vector,
#                               maxSize = 500, minSize = 9, nperm = 10000, positive.enrich.only = FALSE) {
#   require(fgsea)
#   print("finding enrichment at top and bottom of gene ranks, set positive.enrich.only = TRUE if only + enrichment desired")
#   gsea = list()
#   for (i in 1:length(rank.list.celltype)) {
#     gsea[[i]] =
#     fgsea(pathways = pathways,
#     	stats = rank.list.celltype[[i]],
#     	maxSize = maxSize,
#     	minSize = minSize,
#     	nperm = nperm) %>%
#     mutate(celltype = celltypes.vector[i]) %>%
#     arrange(pval)
#     if (positive.enrich.only == TRUE) {
#       gsea[[i]] =
#         gsea[[i]] %>%
#         filter(ES> 0)
#     }
#   }
#   names(gsea) = celltypes.vector
#   return(gsea)
# }
#
#
# # merge GSEA list into dataframe
# RbindGseaResultList = function(gsea_result_list, NES_filter = 0, padj_filter = 0.05){
#
# score = lapply(gsea_result_list, function(x){
#   x = x %>%
#     dplyr::select(pathway, padj, NES, celltype) %>%
#     dplyr::filter( NES > NES_filter )
#   })
# score = do.call(rbind, score) %>%
#   dplyr::filter(padj < padj_filter) %>%
#   dplyr::mutate(n_logp = -log10(padj))
# return(score)
# }




# modules.subset.by.celltype is a list of modules to test with length = celltypes.vector and n modules = unique to subset
# if we use option of only analyzing thresholded modules, u in function above is the length of the number of modules passing threshold.
# AverageSampleModuleZscore = function(average.metacell.list,
#                                      module.list,
#                                      use.module.subset = TRUE,
#                                      modules.subset.by.celltype = modules.test) {
#   mod.scores.celltype = list()
#   if(use.module.subset == TRUE){
#     for (i in 1:length(average.metacell.list)) {
#       mod.scores.celltype[[i]] =
#         calc_avg_module_zscore(module.list = module.list[modules.subset.by.celltype[[i]]],
#                                average.data.frame = average.metacell.list[[i]])
#       }
#     } else {
#     for (i in 1:length(average.metacell.list)) {
#       mod.scores.celltype[[i]] =
#         calc_avg_module_zscore(module.list = module.list[modules.subset.by.celltype[[i]]],
#                                average.data.frame = average.metacell.list[[i]])
#     }
#   }
#   names(mod.scores.celltype) = celltypes
#   return(mod.scores.celltype)
# }




# called sometimes


# usage
# CD4 memory T
# cd4mplot = c("LI.M88.0 leukocyte migration",
#              "LI.M169 mitosis (TF motif CCAATNNSNNNGCG)",
#              "LI.M46 cell division stimulated CD4+ T cells")
# plot.subset.modules(fc.list.entry = 1,
#                     fc.list = fc1,
#                     modules.plot = cd4mplot, width = 3.5,
#                     title.use = "CD4 memory T cells",
#                     filename = "aug2_day1_memcd4_response.png",
#                     path = "btm_analysisV2_jul172019_module_zscores")
#



## LeadingEdge plots
# GetLeadingEdgeGenes = function(gsea.result.list, celltype.index, module.name) {
#
#   celltype = gsea.result.list[[celltype.index]]
#   print("vector of leadingEdge genes for : ")
#   print(unique(celltype$celltype))
#   print(module.name)
#   genes_use =
#     gsea.result.list[[celltype.index]] %>%
#     filter(pathway == module.name) %$%
#     leadingEdge %>%
#     unlist %>%
#     as.character()
#   return(genes_use)
# }
#
#





###################################################
# Modified Seurat version 3 aggretation functions #
###################################################


# This function is a modification on the Seurat AverageExpression. Sum counts for each sample.
# Sum_Counts_Seurat <- function (object, assays = "RNA", features = NULL, return.seurat = FALSE,
#     add.ident = NULL, slot = "counts", use.scale = FALSE, use.counts = FALSE, verbose = TRUE, ...)
# {
#     CheckDots(..., fxns = "CreateSeuratObject")
#     if (use.scale) {
#         .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
#         slot <- "scale.data"
#     }
#     if (use.counts) {
#         .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
#         if (use.scale) {
#             warning("Both 'use.scale' and 'use.counts' were set; using counts",
#                 call. = FALSE, immediate. = TRUE)
#         }
#         slot <- "counts"
#     }
#     fxn.average <- switch(EXPR = slot, data = sum, sum)
#     object.assays <- FilterObjects(object = object, classes.keep = "Assay")
#     assays <- assays %||% object.assays
#     ident.orig <- Idents(object = object)
#     orig.levels <- levels(x = Idents(object = object))
#     ident.new <- c()
#     if (!all(assays %in% object.assays)) {
#         assays <- assays[assays %in% object.assays]
#         if (length(assays) == 0) {
#             stop("None of the requested assays are present in the object")
#         }
#         else {
#             warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
#         }
#     }
#     if (!is.null(x = add.ident)) {
#         new.data <- FetchData(object = object, vars = add.ident)
#         new.ident <- paste(Idents(object)[rownames(x = new.data)],
#             new.data[, 1], sep = "_")
#         Idents(object, cells = rownames(new.data)) <- new.ident
#     }
#     data.return <- list()
#     for (i in 1:length(x = assays)) {
#         data.use <- GetAssayData(object = object, assay = assays[i],
#             slot = slot)
#         features.assay <- features
#         if (length(x = intersect(x = features, y = rownames(x = data.use))) <
#             1) {
#             features.assay <- rownames(x = data.use)
#         }
#         data.all <- data.frame(row.names = features.assay)
#         for (j in levels(x = Idents(object))) {
#             temp.cells <- WhichCells(object = object, idents = j)
#             features.assay <- unique(x = intersect(x = features.assay,
#                 y = rownames(x = data.use)))
#             if (length(x = temp.cells) == 1) {
#                 data.temp <- (data.use[features.assay, temp.cells])
#                 if (slot == "data") {
#                   data.temp <- expm1(x = data.temp)
#                 }
#             }
#             if (length(x = temp.cells) > 1) {
#                 data.temp <- apply(X = data.use[features.assay,
#                   temp.cells, drop = FALSE], MARGIN = 1, FUN = fxn.average)
#             }
#             data.all <- cbind(data.all, data.temp)
#             colnames(x = data.all)[ncol(x = data.all)] <- j
#             if (verbose) {
#                 message(paste("Finished averaging", assays[i],
#                   "for cluster", j))
#             }
#             if (i == 1) {
#                 ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
#             }
#         }
#         names(x = ident.new) <- levels(x = Idents(object))
#         data.return[[i]] <- data.all
#         names(x = data.return)[i] <- assays[[i]]
#     }
#     if (return.seurat) {
#         toRet <- CreateSeuratObject(counts = data.return[[1]],
#             project = "Sum", assay = names(x = data.return)[1],
#             ...)
#         if (length(x = data.return) > 1) {
#             for (i in 2:length(x = data.return)) {
#                 toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
#             }
#         }
#         if (DefaultAssay(object = object) %in% names(x = data.return)) {
#             DefaultAssay(object = toRet) <- DefaultAssay(object = object)
#         }
#         Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
#         Idents(object = toRet) <- factor(x = Idents(object = toRet),
#             levels = as.character(x = orig.levels), ordered = TRUE)
#         toRet <- NormalizeData(object = toRet, verbose = verbose)
#         toRet <- ScaleData(object = toRet, verbose = verbose)
#         return(toRet)
#     }
#     else {
#         return(data.return)
#     }
# }

# environment(Sum_Counts_Seurat) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions

# modified version of Seurat's AverageExpression function to get the average of (in case of Flu analysis)the scran normalized counts.
# This makes an average "meta cell"dataset. where each count contributes equally to the aggregated library.
# designed to aggregate on a subject level and in the DE pipeline usage is to make a average list of objects per subject by cell type.
# example usage

# source("AverageExpresson2.r")
# al = list()
# for (i in 1:length(celltypes)) {
#   # get cells in each cell tpe
#   cells =
#     h1 %>%
#     SetAllIdent(id = "cell_type_ADT_0.3") %>%
#     WhichCells(ident = celltypes[i])
#   al[[i]] =
#     SubsetData(h1, cells.use = cells) %>%
#     AverageExpression2(add.ident = "sample", use.raw = F, use.scale = F, return.seurat = F)
# }
# names(al) = celltypes
# saveRDS(al, file = "../cell_subset_average_list.rds")

# above sample was e.g. 200_d1 which is subject 200 day 1 cells by celltype.
# only average the gene expression data not the other assay data.
# change to taking the mean ofte lg normalized counts not exp - 1 of the normalized counts .
# AverageExpression3 <- function (object, assays = NULL, features = NULL, return.seurat = FALSE,
#     add.ident = NULL, slot = "data", use.scale = FALSE, use.counts = FALSE,
#     verbose = TRUE, ...)
# {
#     CheckDots(..., fxns = "CreateSeuratObject")
#     if (use.scale) {
#         .Deprecated(msg = "'use.scale' is a deprecated argument, please use the 'slot' argument instead")
#         slot <- "scale.data"
#     }
#     if (use.counts) {
#         .Deprecated(msg = "'use.counts' is a deprecated argument, please use the 'slot' argument instead")
#         if (use.scale) {
#             warning("Both 'use.scale' and 'use.counts' were set; using counts",
#                 call. = FALSE, immediate. = TRUE)
#         }
#         slot <- "counts"
#     }
#     fxn.average <- switch(EXPR = slot, data = mean, mean)
#     object.assays <- FilterObjects(object = object, classes.keep = "Assay")
#     assays <- assays %||% object.assays
#     ident.orig <- Idents(object = object)
#     orig.levels <- levels(x = Idents(object = object))
#     ident.new <- c()
#     if (!all(assays %in% object.assays)) {
#         assays <- assays[assays %in% object.assays]
#         if (length(assays) == 0) {
#             stop("None of the requested assays are present in the object")
#         }
#         else {
#             warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
#         }
#     }
#     if (!is.null(x = add.ident)) {
#         new.data <- FetchData(object = object, vars = add.ident)
#         new.ident <- paste(Idents(object)[rownames(x = new.data)],
#             new.data[, 1], sep = "_")
#         Idents(object, cells = rownames(new.data)) <- new.ident
#     }
#     data.return <- list()
#     for (i in 1:length(x = assays)) {
#         data.use <- GetAssayData(object = object, assay = assays[i],
#             slot = slot)
#         features.assay <- features
#         if (length(x = intersect(x = features, y = rownames(x = data.use))) <
#             1) {
#             features.assay <- rownames(x = data.use)
#         }
#         data.all <- data.frame(row.names = features.assay)
#         for (j in levels(x = Idents(object))) {
#             temp.cells <- WhichCells(object = object, idents = j)
#             features.assay <- unique(x = intersect(x = features.assay,
#                 y = rownames(x = data.use)))
#             if (length(x = temp.cells) == 1) {
#                 data.temp <- (data.use[features.assay, temp.cells])
#                 if (slot == "data") {
#                   data.temp <- expm1(x = data.temp)
#                 }
#             }
#             if (length(x = temp.cells) > 1) {
#                 data.temp <- apply(X = data.use[features.assay,
#                   temp.cells, drop = FALSE], MARGIN = 1, FUN = fxn.average)
#             }
#             data.all <- cbind(data.all, data.temp)
#             colnames(x = data.all)[ncol(x = data.all)] <- j
#             if (verbose) {
#                 message(paste("Finished averaging", assays[i],
#                   "for cluster", j))
#             }
#             if (i == 1) {
#                 ident.new <- c(ident.new, as.character(x = ident.orig[temp.cells[1]]))
#             }
#         }
#         names(x = ident.new) <- levels(x = Idents(object))
#         data.return[[i]] <- data.all
#         names(x = data.return)[i] <- assays[[i]]
#     }
#     if (return.seurat) {
#         toRet <- CreateSeuratObject(counts = data.return[[1]],
#             project = "Average", assay = names(x = data.return)[1],
#             ...)
#         if (length(x = data.return) > 1) {
#             for (i in 2:length(x = data.return)) {
#                 toRet[[names(x = data.return)[i]]] <- CreateAssayObject(counts = data.return[[i]])
#             }
#         }
#         if (DefaultAssay(object = object) %in% names(x = data.return)) {
#             DefaultAssay(object = toRet) <- DefaultAssay(object = object)
#         }
#         Idents(toRet, cells = colnames(x = toRet)) <- ident.new[colnames(x = toRet)]
#         Idents(object = toRet) <- factor(x = Idents(object = toRet),
#             levels = as.character(x = orig.levels), ordered = TRUE)
#         toRet <- NormalizeData(object = toRet, verbose = verbose)
#         toRet <- ScaleData(object = toRet, verbose = verbose)
#         return(toRet)
#     }
#     else {
#         return(data.return)
#     }
# }

# environment(AverageExpression3) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions


# # get vector of unique celltypes
# GetCelltypes = function(SeuratObject, celltype) {
#   md = SeuratObject@meta.data
#   x = md[[celltype]] %>% unique()
#   return(x)
# }




# create a list of summarized data aggregated by summing raw counts by sample (subject1_time1) and celltype
# slot = counts gets raw data.
# MakePseudobulkList = function(seurat.object, celltype_column, sample_column, vector.of.celltypes) {
#   sl = list()
#   Idents(seurat.object) <- celltype_column
#   for (i in 1:length(vector.of.celltypes)) {
#       seurat.object.temp <- subset(seurat.object, idents = vector.of.celltypes[i])
#       sl[[i]] = Sum_Counts_Seurat(seurat.object.temp, add.ident = sample_column,
#       	assays = "RNA",
#       	features = NULL,
#       	return.seurat = FALSE,
#       	slot = "counts",
#       	verbose = TRUE)
#   }
#   names(sl) = vector.of.celltypes
#   return(sl)
# }


## Same as above but it calls AverageExpression3
# MakeMetaCellList = function(seurat.object, celltype_column, sample_column, vector.of.celltypes) {
#   sl = list()
#   Idents(seurat.object) <- celltype_column
#   for (i in 1:length(vector.of.celltypes)) {
#       seurat.object.temp <- subset(seurat.object, idents = vector.of.celltypes[i])
#       sl[[i]] = AverageExpression3(seurat.object.temp,
#       	add.ident = sample_column,
#       	assays = "RNA",
#       	features = NULL,
#       	return.seurat = FALSE,
#       	slot = "data",
#       	verbose = TRUE)
#   }
#   names(sl) = vector.of.celltypes
#   return(sl)
# }



# # GEt a metadata table
# MakeMetaTableFromSeurat = function(seurat_object, sample_column, variable_column,
#                                    aggregate.data.list, celltypes.vector = celltypes){
#   meta = seurat_object@meta.data %>% select(sample_column, variable_column)
#   meta = table(meta[[sample_column]], meta[[variable_column]]) %>%
#     as.data.frame.table() %>%
#     filter(Freq > 0 ) %>%
#     select(-Freq) %>%
#     column_to_rownames("Var1")
#   names(meta) = variable_column
#
#   # match up the order of the metadata table to the aggregate data frames.
#   matchcols =
#     colnames(aggregate.data.list[[1]]$RNA) %>% as.data.frame()
#   names(matchcols) = "sample"
#   matchcols =
#     matchcols %>%
#     mutate_if(is.factor, as.character()) %>%
#     #mutate(sample_name = str_replace(sample, pattern = paste0(celltypes[1],"_"), replacement = ""))
#     mutate(sample_name = str_replace(sample, pattern = paste0(celltypes.vector[1],"_"), replacement = ""))
#
#   meta = meta %>% rownames_to_column("sample")
#   meta = meta[match(matchcols$sample_name, meta$sample), ]
#   return(meta)
# }


### Avg module score functions.

# run on each cell subset single cell level data:
# SeuratListByCelltype = function(seurat.object, celltype.ident, return.matrix = T, celltypes.vector){
#   seurat.object = Idents(seurat.object) <- celltype.ident
#   sl = list()
#   for (i in 1:length(celltypes.vector)) {
#     sl[[i]] =  subset(seurat.object, idents = celltypes.vector[i])
#     if (return.matrix == TRUE) {
#       sl[[i]] = sl[[i]]@data %>% as.matrix
#     } else {
#       print(class(sl[[i]]))
#     }
#   }
#   return(sl)
# }
# usage
## these are the celltypes defined from the get.subject.celltype.table function
# celltypes = unique(h1@meta.data$celltype_V2)
# celltypes = celltypes[-c(2,10,15,16,17,18,19,21,22,23)]
#sl = SeuratListByCelltype(seurat.object = h1, celltype.ident = "celltype_V2", return.matrix = T, celltypes.vector = celltypes)


######## AUCell methods
# returns the list of cell by gene matrices by (protein defined) cell-type above filtered:
# each celltype retains genes expressed in at least "percent.express.threshold" percent of cells, default 1%
# FilterGenesByCelltype = function(list.by.celltype, percent.express.threshold = 0.01, celltypes.vector){
#   filtered.matrix = list()
#   for (i in 1:length(list.by.celltype)) {
#     # percent of cells that the gene is expressed in
#     genes_keep = apply(list.by.celltype[[i]], 1, function(x) { length(which(x > 0)) / length(x) })
#     # subset by a percent threshold (here 1% )
#     genes_use = genes_keep[genes_keep > percent.express.threshold] %>% names
#     print(length(genes_use))
#     filtered.matrix[[i]] = list.by.celltype[[i]][genes_use, ]
#   }
#   names(filtered.matrix) = celltypes.vector
#   return(filtered.matrix)
# }

# get modules that have at least 20% of their genes in the filtered data
# ModuleAUCRankFilter = function(filtered.matrix, module.list, celltypes.vector = celltypes,
#                                ncores = 6, save.ranks.path = dir.create("auc_ranks_output")){
#   require(AUCell)
#   cellrank = cellauc = modules.test = list()
#   if (!exists(x = "auc_ranks_output")) {
#     dir.create(path = "auc_filter_output")
#   }
#   for (i in 1:length(filtered.matrix)) {
#     print("calculating AUC stats for")
#     print(celltypes.vector[i])
#     sink(paste0("auc_filter_output/",celltypes.vector[i],"ranks.txt"))
#     cellrank[[i]] = AUCell_buildRankings(exprMat = filtered.matrix[[i]], plotStats = F, nCores = 8, verbose = T)
#     sink()
#     sink(paste0("auc_filter_output/",celltypes.vector[i],"AUC.txt"))
#     cellauc[[i]]= AUCell_calcAUC(geneSets =  btm, rankings = cellrank[[i]],
#                                  aucMaxRank = nrow(cellrank[[i]])*0.20, nCores = 8, normAUC = T)
#     sink()
#     # modules with 20% representation in the filtered dataset
#     modules.test[[i]] = cellauc[[i]]@NAMES
#     saveRDS(cellrank[[i]], file = paste0(save.ranks.path,celltypes.vector[i],".rds"))
#     saveRDS(cellauc[[i]], file = paste0(save.ranks.path,celltypes.vector[i],".rds"))
#   }
#   print("AUC statistics for each cell type are in dir: auc_filter_output")
#   names(modules.test) =  celltypes.vector
#   return(modules.test)
# }


### depends AverageExpression3


