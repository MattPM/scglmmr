#' #' dream.mixed.model - run dream mixed model
#' #' @param dge_lists list of dgelists created with NormalizePseudobulk
#' #' @param apriori_contrasts one of TRUE or FALSE, whether to fit a priori contrasts
#'
#' #' @param contrast_matrix contrast matrix created with make.contrasts
#' #' @param design_matrix design matrix created with BulkDesignMatrix
#' #' @param plotsavepath a path to save created plot of contrasts
#' #' @param ncores number of cores for doParallel
#' #' @param cell_metadata metadata cells-rows variables-columns i.e. ColData or seuratmeta.data
#' #' @param fixed_effects a vector of covariates that are columns in metadata e.g. a vector of  'age' 'gender'
#' #' @param sample_column quoted character e.g. "sample" the subject level sample variable should have multiple timepoints subjectID_timepoint i.e. s1_0, s1_1
#' #' @param lme4_formula symbolic model formula the default is '~ 0 + cohort_timepoint + (1|sampleid)'
#' #'
#' #' @return list of model fits indexed by celltype
#' #' @importFrom variancePartition voomWithDreamWeights plotContrasts dream
#' #' @importFrom limma voom
#' #' @importFrom parallel makeCluster
#' #' @importFrom BiocParallel register SnowParam
#' #' @importFrom rlang sym
#' #' @importFrom dplyr group_by summarize_each
#' #' @export
#' #'
#' #' @examples
#' #'\dontrun{
#' #'c_mat = makeContrasts(
#' #'foldchange_difference = (group_timepoint1_1 - group_timepoint1_0) - (group_timepoint0_1 - group_timepoint0_0),
#' #'time1_foldchange = (group_timepoint1_1 + group_timepoint0_1) / 2  - (group_timepoint1_0 + group_timepoint0_0) / 2,
#' #'baseline_groups = (group_timepoint1_0 - group_timepoint0_0),
#' #'levels = colnames(designmat)
#' #')
#' #'# fit mixed model for the multi timepoint contrasts
#' #'fit = scglmmr::dreamMixedModel(dge_lists = dge,
#' #'apriori_contrasts = TRUE,
#' #'sample_column = 'sample',
#' #'cell_metadata = meta,
#' #'contrast_matrix = c_mat,
#' #'design_matrix = designmat,
#' #'lme4_formula =  '~ 0 + age + gender + cohort_timepoint + (1|sampleid)',
#' #'fixed_effects = c('age', 'gender', 'cohort_timepoint'),
#' #'plotsavepath = figpath,
#' #'ncores = 4)
#' #' }
#' #' # run dream mixed model
#' # dream method Hoffman et.al. 2020  https://doi.org/10.1093/bioinformatics/btaa687
#' dreamMixedModel = function(dge_lists, apriori_contrasts = FALSE,
#'                            sample_column,
#'                            sample.metadata,
#'                            contrast_matrix = NULL,
#'                            design_matrix,
#'                            fixed_effects,
#'                            cell_metadata,
#'                            lme4_formula = '~ 0 + cohort_timepoint + (1|sampleid)',
#'                            ncores = 4, ...) {
#'
#'   BiocParallel::register(BiocParallel::SnowParam(workers = ncores))
#'   pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
#'   # cl = parallel::makeCluster(ncores)
#'   # doParallel::registerDoParallel(cl = ncores)
#'
#'   if(isFALSE('subjectID' %in% colnames(cell_metadata))){
#'     stop("metadata must have a column 'sampleid' to fit a varying intercept model with subjectIDs")
#'   }
#'
#'   # combine vector of covariates in model
#'   vars_all = c(sample_column, fixed_effects, 'sampleid')
#'
#'   # define sample_column argument as data-variable to make non environment variable
#'   gvar = rlang::sym(sample_column)
#'   # collapse cells down to data.variables in model formula {{}} for non standard evaluation of data.variable
#'   model_md = cell_metadata[ ,vars_all] %>%
#'     dplyr::group_by({{gvar}}) %>%
#'     dplyr::summarise_each(list(~unique(.)))
#'
#'   # convert tibble to dataframe to use in dream fit
#'   model_md = base::as.data.frame(model_md)
#'   rownames(model_md) = model_md[[sample_column]]
#'   print("dream data argument for model "); print(model_md)
#'   # check row index of model metadata matches dgelist columns and show symbolic formula
#'   stopifnot(isTRUE(all.equal(target = colnames(dge_lists[[1]]), current = rownames(model_md))))
#'
#'   print('model specified (change with argument to lme4_formula) '); print(lme4_formula)
#'   v1 = lapply(dge_lists, function(x){
#'     design = NULL
#'     variancePartition::voomWithDreamWeights(counts = x, data = model_md, formula = lme4_formula,
#'                                             normalize.method = "none", save.plot = TRUE, plot = TRUE)
#'   })
#'   if(isTRUE(apriori_contrasts)) {
#'     # fit mixed model
#'     fit1 = lapply(v1, function(x){
#'       variancePartition::dream(...,
#'                                exprObj = x,
#'                                formula = lme4_formula,
#'                                data = model_md,
#'                                L = as.matrix(contrast_matrix),
#'                                BPPARAM = pparam)
#'     })
#'   } else {
#'     # fit model without custom contrasts
#'     fit1 = lapply(v1, function(x){
#'       variancePartition::dream(exprObj = x, formula = lme4_formula, data = model_md)
#'     })
#'   }
#'   return(fit1)
#' }
#'
#'
#'
#' AggregateCellMetadata = function(cell.metadata, variables, sample.column){
#'
#'   if (!sample.column %in% colnames(cell.metadata)) {
#'     stop('sample column was not found in column names of cell metadata')
#'   }
#'   if (!variables %in% colnames(cell.metadata)) {
#'     stop('some `variables` were not found in column names of cell metadata')
#'   }
#'
#'   # combine vector of covariates in model
#'   vars_all = c(variables, sample.column)
#'
#'   # define sample_column argument as data-variable to make non environment variable
#'   gvar = rlang::sym(sample.column)
#'   # collapse cells down to data.variables in model formula
#'   model_md = cell_metadata[ ,vars_all] %>%
#'     dplyr::group_by({{gvar}}) %>%
#'     dplyr::summarise_each(list(~unique(.)))
#'
#'   # convert tibble to dataframe to use in dream fit
#'   model_md = base::as.data.frame(model_md)
#'   rownames(model_md) = model_md[[sample.column]]
#'   print("aggregated metadata summary "); print(head(model_md))
#' }
#'
#'
#'
#' # example of using a custom contrast:
#' # re-level time.group into ordered combined factor
#' samplemd$time.group = factor(
#'   samplemd$time.group,
#'   levels = c('d0_group1', 'd1_group1', 'd0_group2', 'd1_group2')
#'   )
#'
#' # specify random intercept model
#' f1 <- ~ 0 + time.group + gender + scaledage + (1|subjectid)
#'
#' # specify contrast matrix to test the fold change difference
#' # based on levels of time.group this should be cmat = c(-1, 1, 1, -1, 0, 0)
#' # the age and gender other main effects in the model should all be 0 in L
#' L2 = makeContrastsDream(
#'   formula = f1,
#'   data = samplemd,
#'   contrasts = c(
#'     baseline = "time.groupd0_group1- time.groupd0_group2",
#'     delta = "(time.groupd1_group1 - time.groupd0_group1) - (time.groupd1_group2 - time.groupd0_group2)",
#'     treatment = "( time.groupd1_group1 + time.groupd1_group2 ) / 2 - ( time.groupd0_group1 + time.groupd0_group2 ) / 2 "
#'   )
#' )
#'
#'
#' dream::plotContrasts(L2) +
#'   ggsave(filename = paste0(figpath,'contrastmodel.pdf'),
#'          width = 7, height = 4)
#'
#' # fit model on each subset
#' # init store
#' fit1 = v1 = list()
#' for (i in 1:length(pb)) {
#'
#'   # init data
#'   meta = samplemd
#'   form = f1
#'   contrast_matrix = L2
#'   counts = pb[[i]]
#'
#'   # dge list
#'   d = edgeR::DGEList(counts = counts, samples = meta)
#'
#'   # filter cell type specific lowly expressed genes and calc norm factors
#'   gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = as.factor(d$samples$time.group))
#'   print(names(pb)[i]);print(table(gtable))
#'   d = d[gtable, keep.lib.sizes=FALSE]
#'   d = edgeR::calcNormFactors(object = d)
#'
#'   # get voom observation level weights
#'   v = voomWithDreamWeights(counts = d,
#'                            formula = form,
#'                            data = meta,
#'                            BPPARAM = pparam,
#'                            plot = TRUE, save.plot = TRUE)
#'   # fit contrast mixed model
#'   fitmm = dream(exprObj = v,
#'                 formula = form,
#'                 data = meta,
#'                 L = contrast_matrix,
#'                 BPPARAM = pparam,
#'                 useWeights = TRUE, REML = TRUE)
#'   # save results
#'   v1[[i]] = v
#'   fit1[[i]] = fitmm
#' }
#' names(v1) = names(fit1) = names(pb)
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
