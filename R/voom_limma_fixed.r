# scglmmr pseudobulk differential expression pipeline
# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

#' RunVoomLimma
#'
#' @param dgelists a ist of DGElist created in NormalizePseudobulk
#' @param design_matrix design matrix
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
    limma::voom(counts = x,
                design = design_matrix,
                normalize.method = "none",
                save.plot = FALSE,
                plot = TRUE)
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

################# AJM function

# RunVoomLimma_covar - Note:
# removed celltypes vector argument to match updated function above; argument was not needed.
# same as above [RunLimmaVoom] with option to adjust for covariate within the mRNA


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
  return(eb_c_fit)
}
