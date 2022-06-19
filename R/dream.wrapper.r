# scglmmr pseudobulk differential expression pipeline
# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com



#' FitDream - run mixed effects model on aggregated (summed) data using the method 'dream'
#' by Hoffman et. al. Bioinformatics (2021) doi.org/10.1093/bioinformatics/btaa687. Fits
#' mixed model using lme4 with REML and voom weights.
#' @param dge.lists list of DGEList objects indexed by cell types -- the object returned by `scglmmr::Normalize`
#' @param sample.metadata metadata, for example, object returned by AggregateCellMetadata
#' @param pparam number of cores for biocparallel. Set with BiocParallel::register(BiocParallel::SnowParam(4)); pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE). use the desired number of cores
#' @param sample_column quoted character e.g. "sample" the subject level sample variable should have multiple timepoints subjectID_timepoint i.e. s1_0, s1_1
#' @param lme4.formula symbolic model formula for model to be fit, for example,
#' '~ 0 + group.timepoint + age + sex + (1|SubjectID)'. covariates must be in sample.metadata
#'
#' @return list of model fits indexed by celltype
#' @importFrom variancePartition voomWithDreamWeights dream eBayes
#' @importFrom BiocParallel register SnowParam
#'
#' @export
#'
#' @examples
#'\dontrun{
#'# make contrast matrix
#' L2 = makeContrastsDream(
#'   formula = f1,
#'   data = metadata,
#'   contrasts = c(
#'     baseline = "Group.time1_0 - Group.time0_0",
#'     treatment_delta = "( Group.time1_1 - Group.time1_0 ) - ( Group.time0_1 - Group.time0_0 )",
#'     treatment = "( Group.time1_1 + Group.time0_1 ) / 2 - ( Group.time1_0 + Group.time0_0 ) / 2 "
#'   )
#' )
#' f1 = '0 + group.time + age + sex + (1|SubjectID)'
#'
#' }
#'
#' fits = FitDream(pb.list = pb, sample.metadata = metadata, lme4.formula = f1, dream.contrast.matrix = L2, ncores = 4)

FitDream = function(pb.list,
                    sample.metadata,
                    lme4.formula,
                    dream.contrast.matrix = NULL,
                    returnvoom = FALSE,
                    ncores = 4, ...){


  print(' Fitting models with dream method ')
  print(' If using this model cite Hoffman et. al. Bioinformatics (2021) doi.org/10.1093/bioinformatics/btaa687')

  # checks
  if (is.null(lme4.formula)) {
    stop('specify lme4.formula e.g. f <- ~ 0 + time.group + sex + age + (1|subjectid)')
  }
  if(is.null(pparam)){
    stop('specify number of cores for multi threading. Prior to running FitDream, use BiocParallel::register(BiocParallel::SnowParam(4)); pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE). example - 4 ores, use the desired number of cores')
  }

  # init data
  meta = sample.metadata
  form = lme4.formula
  pb.list = dge_lists

  # init store
  res = list()
  for (i in 1:length(pb.list)) {
    d = pb.list[[i]]


    print(
      paste0("fitting mixed models for ",
             print(names(pb.list)[i])," ",i, " of ", length(pb.list), "subsets")
          )

    # calculate voom observation level weights
    v = voomWithDreamWeights(
      counts = d,
      formula = form,
      data = meta,
      BPPARAM = pparam,
      plot = TRUE,
      save.plot = TRUE
    )

    # fit models
    if (!is.null(contrast_matrix)) {
      fitmm = dream(...,
                    exprObj = v,
                    formula = form,
                    data = meta,
                    L = contrast_matrix,
                    BPPARAM = pparam,
                    useWeights = TRUE,
                    REML = TRUE
                    )
      fitmm = variancePartition::eBayes(fitmm)
    } else{
      fitmm = dream(...,
                    exprObj = v,
                    formula = form,
                    data = meta,
                    BPPARAM = pparam,
                    useWeights = TRUE,
                    REML = TRUE
                    )
      fitmm = variancePartition::eBayes(fitmm)
    }

    # return
    if(isTRUE(returnvoom)){
      res[[i]] = list('fit' = fitmm, 'voom' = v)
    }else{
      res[[i]] = fitmm
    }
  } # celltype i
  names(res) = names(pb.list)
}

