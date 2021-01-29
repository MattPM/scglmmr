#' WeightedCellModuleScore - calculate the average or weighted aveage +/- scaling return dataframe of cells by average module score
#'
#' @param seurat_object - a Seurat object - if using version 2 specify `version` = "2", if no seurat object specified, gene_matrix and cell_metadata must be specified
#' @param gene_matrix  - normalized genes (rows) by cells (columns)
#' @param module_list - names list of gene modules - each element a vector of gene names.
#' @param threshold - at least this fraction of genes in the signature must be < 0 across all cells or else a score is not calculated
#' @param cellwise_scaling - T/F scale across cells ?
#' @param return_weighted - T/F weight the averae by multiplying by gene representation?
#' @param Seurat_version - if seurat_object is not NULL, and seurat object version < 3 input "2"
#'
#' @return
#' @importFrom Matrix rowSums colMeans
#' @import Seurat
#' @export
#'
#' @examples
WeightedCellModuleScore = function(seurat_object= NULL, gene_matrix = NULL, module_list,
                                   threshold = 0.1, cellwise_scaling = FALSE, return_weighted = TRUE, Seurat_version = NULL){

  # may remove this section to not require import of Seurat but keeping for now for backwards compatibility with FSC
  # require(Seurat)

  if(!is.null(seurat_object)) {
    message("extracting gene expression from seurat object")
    if(Seurat_version == "2" & !is.null(seurat_object@data)) {
      mtx = seurat_object@data
    } else if(!is.null(seurat_object@assays$RNA@data)) {
      mtx = seurat_object@assays$RNA@data
    } else {
    mtx = gene_matrix
    }
  }
  ###################################


  # calculate score for each module
  for (i in 1:length(module_list)) {
    # init
    signature = module_list[[i]]
    signaturename = names(module_list[i])
    # calc weights
    n_modules = length(module_list)
    print(paste0("calculating module score for module ",i, " of ",n_modules))
    # get expressed genes
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
        # calculate av base score to be modified by return args
        mod_avg = Matrix::colMeans(mtx[gene_universe, ])
      }
    # return options
    # option to return score scaled across all cells / can change to return non weighted score
    if (isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell z scored of cells in object
      score_return = scaled_mod_weight_avg = scale(mod_avg * frac_nonzero)
    }
    if (!isTRUE(cellwise_scaling) & isTRUE(return_weighted)) {
      # average weighted by non-zero gene representation in cell
      score_return = weighted_average = mod_avg * frac_nonzero
    }
    if (isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average in cell z scored of cells in object
      score_return = scaled_mod_avg =  scale(mod_avg)
    }
    if (!isTRUE(cellwise_scaling) & !isTRUE(return_weighted)) {
      # average exprs of genes in cell
      score_return =  mod_avg
    }
    # format returned score as dataframe columns = modules, rownames = barcodes
    score_return = as.data.frame(score_return)
    names(score_return) = names_vec[i]
    score_keep[[i]] = score_return
  }
  score_df = do.call(cbind, score_keep)
  return(score_df)
  # can add this back to Seurat Object or SCE object as metadata
}
