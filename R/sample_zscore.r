# scglmmr pseudobulk differential expression pipeline
# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com


#' calc_avg_module_zscore - calculate average module z score of list of modules on a PseudobulkList
#' This is equivalent to the average z score method used in in Kotliarov et. al. Nature Med 2020
#' zscore is calculated across both genes and samples
#' it is adopted below to run on 'pseudobulk lists' (average "averagemetacell.list" or pseudobulk list
#' created by PseudobulkList) this small wrapper is called by the AverageSampleModuleZscore.
#' calculate signature score for each cell type, BTM, Subject
#' function input = named list of modules, dataframe with subject as rows genes as columns
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


#' AverageSampleModuleZscore apply the function calc_avg_module_zscore to a pseudobulklist.
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
#' # av  = object returned by PseudobulkList
#' av.zscore = AverageSampleModuleZscore(
#'  average.metacell.list = av,
#'  module.list = btm,
#'  use.module.subset = FALSE,
#'  )
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
