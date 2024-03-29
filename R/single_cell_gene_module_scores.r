# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

#' WeightedCellModuleScore - calculate the average or weighted aveage +/- scaling return dataframe of cells by average module score
#'
#' @param gene_matrix  - normalized genes (rows) by cells (columns)
#' @param module_list - names list of gene modules - each element a vector of gene names.
#' @param threshold - at least this fraction of genes in the signature must be < 0 across all cells or else a score is not calculated
#' @param cellwise_scaling - T/F scale across cells ?
#' @param return_weighted - T/F weight the averae by multiplying by gene representation?
#'
#' @return dataframe of cells barcodes (rownames) by gene module scores (columns)
#' @importFrom Matrix rowSums colMeans
#' @export
#'
#' @examples
#'\dontrun{
#
# # load data from single cell data object
# Seurat = readRDS("my_seurat_object.rds")
#
# # add cellwise module score for each signature
# mod_scores = WeightedCellModuleScore(seurat_object = Seurat,
#                                      module_list = btm,
#                                      threshold = 0.1,
#                                      return_weighted = FALSE, cellwise_scaling = FALSE,
#                                      Seurat_version = "2")
# Seurat = AddMetaData(Seurat,metadata = mod_scores)
# module_n = names(sig_test)
#
# # set up module data frame
# module_df = Seurat@meta.data %>% select(barcode_check, celltype_joint, module_n)
#
# # format metadata as factors group_id is order leveled for contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
# md = Seurat@meta.data %>%
#   mutate(group_id = factor(treat_time,  levels = c('pre_low', 'post_low', 'pre_high', 'post_high'))) %>%
#   mutate(sampleid = factor(sampleid)) %>%
#   select(barcode_check, celltype_joint, sampleid,  age, group_id)
#
# # Fit mixed model
# plot_savepath = paste0(my_figure_save_path, "/marginalmeans/"); dir.create(plot_savepath)
#
# # specify any random intercept model e.g.
# f1 = 'modulescore ~ age + group_id + (1|sampleid)'
#
# # fit sc mod mixed model on ewighted module scores.
# mm_res = SCGroupContrastGLMM(module_data_frame = module_df,
#                              celltype_column = 'celltype',
#                              metadata = md,
#                              fixed_effects = NULL,
#                              lmer_formula = f1,
#                              plotdatqc = TRUE,
#                              figpath = 'your/file/path')
#'}
WeightedCellModuleScore = function(gene_matrix = NULL,module_list,threshold = 0,
                                   cellwise_scaling = FALSE,return_weighted = FALSE){


  mtx = gene_matrix
  score_keep = list()
  for (i in 1:length(module_list)) {
    # init storage
    signature = module_list[[i]]
    # calc weights
    signaturename = names(module_list[i])
    n_modules = length(module_list)
    print(paste0("calculating module score for module ",i, " of ",n_modules))
    gene_universe = intersect(signature, rownames(mtx))

    # extract number of genes with any expression across cells
    rs = Matrix::rowSums(mtx[gene_universe, ])
    frac_nonzero = length(signature[rs>0]) / length(signature)

    # do not calculate module score for modules with expression less than the param `threshold`
    if( frac_nonzero < threshold ) {
      print(paste0(" fraction of genes in signature with non zero expression = ",  frac_nonzero ,
                   " for " , signaturename, " not calculating score; set threshold to 0 to override",
                   " current threshold = ",  threshold))
      i = i + 1
      } else {
        # calculate the fraction of genes in the signature with nonzero expression in each cell
        gene_rep = apply(mtx[gene_universe, ], 2, FUN = function(x){ length(x[x>0]) / length(signature) } )
        # calculate av score that can be modified by return args
        mod_avg = Matrix::colMeans(mtx[gene_universe, ])
      }
    # return opts
    if (isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell z scored of cells in object
      score_return = scaled_mod_weight_avg = scale(mod_avg * frac_nonzero)
    }
    if (!isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell
      score_return = weighted_average = mod_avg * frac_nonzero
    }
    if (isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average in cell z scored of cells
      score_return = scaled_mod_avg = scale(mod_avg)
    }
    if (!isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average exprs of genes in cell
      score_return =  mod_avg
    }
    score_return = as.data.frame(score_return)
    names(score_return) = signaturename
    score_keep[[i]] = score_return
  }
  score_df = do.call(cbind, score_keep)
  return(score_df)
}
