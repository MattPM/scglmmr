#' SCMixedPoisson
#'
#' @param gene_data matrix of raw single cell UMI counts genes as COLUMNS and cells as ROWS; i.e. a transpose t() of Bioconductor or Seurat raw.data slots if genes_test is specified, columns must contain all genes from genes_test i.e. umi_matrix = t(seurat_raw_counts[unique(unlist((genes_test)), ])
#' @param lme4metadata metadata for each cell and cell barcodes as
#' @param model_formula defaults to 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)' the variable 'timepoint' can be altered to be any treatment or perturbation effect being tested. Other variables should be kept the same.
#' @param reduced_model_formula identical formula as `full_model_formula` without the treatment or perturbation variable term. Defaults to 'gene ~ offset(log(nUMI)) + (1|subjectid)'
#' @param celltype_genes_test A R list indexed by cell type of the subset of genes to test for each cell type.
#'
#' @return a results matrix for all genes across all cell types tested.
#' @export
#'
#' @importFrom lme4 glmer
#'
#' @examples
#' results = SCMixedPoisson(gene_data = gene_data,
#'   lme4metadata = meta_data,
#'   model_formula = 'gene ~ offset(log(nUMI)) + timepoint + (1|sampleid)',
#'   celltype_genes_test = celltype_indexed_gene_vector
#' )
#'
#'
#'
SCMixedPoisson = function(gene_data,
                          lme4metadata,
                          model_formula = 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)',
                          # reduced_model_formula = 'gene ~ offset(log(nUMI)) + (1|subjectid)',
                          celltype_genes_test = NULL){

  ############# test function
  # model_formula = 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)'
  # reduced_model_formula = 'gene ~ offset(log(nUMI)) + (1|subjectid)'
  # gene_data = gene_data
  # lme4metadata = md
  # celltype_genes_test = genes_test
  #############

  # function start
  f1 = model_formula

  ## Run pre model fit checks
  ctmd = split(lme4metadata,f = lme4metadata$celltype)
  celltype_test =
    do.call(rbind,
            lapply(ctmd, function(x){
              any(as.matrix(table(paste(x$sampleid, x$timepoint))) <=1)
            })
    )
  test_out = celltype_test[ ,1][celltype_test[ ,1] == TRUE]
  if (!length(test_out) == 0) {
    stop(paste0('some subjects have 1 or 0 cells for celltype',
                names(celltype_test[ ,1])[celltype_test[ ,1] == TRUE])
    )
  }
  if (!c("celltype") %in% colnames(lme4metadata)) {
    stop("in lme4metadata, rename the cell type column 'celltype' and the random effect (1|subjectid)")
  }
  # Run checks for offset term that should be specified in Poisson model
  if ( isFALSE(grepl("offset(log(nUMI))", f1, fixed = TRUE))) {
    stop("include offset(log(numinUMI)) in f1 and f2")
  }
  if ( isFALSE(grepl("(1|subjectid)", f1, fixed = TRUE))) {
    stop("include (1|subjectid) as a random effect in f1 and f2")
  }
  if ( isFALSE(grepl("offset(log(nUMI))", f2, fixed = TRUE))) {
    stop("include offset(log(numinUMI)) in f1 and f2")
  }
  # Run checks for minimal required structure of f1 and f2 lme4 formulae
  if ( isFALSE(grepl("gene" , f1, fixed = TRUE))) {
    stop(" 'gene' must be specified in LHS of f1 e.g. 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)'")
  }
  if ( isFALSE(grepl("gene" , f2, fixed = TRUE))) {
    stop("'gene' must be specified in LHS of f2 e.g. 'gene ~ offset(log(nUMI)) + (1|subjectid)'")
  }
  # print(paste0("testing poisson GLMM of ", f1, " vs reduced model ", f2, ' with log link'))

  # init storage
  res_celltype = res_list = list()

  # define gene names to test
  gene_names = colnames(gene_data)

  # define celltypes
  cts = as.character(unique(lme4metadata$celltype))

  # make sure celltypes ordered by subset of genes test list by cell type
  if(!is.null(celltype_genes_test)){
    stopifnot(all.equal(names(celltype_genes_test), cts))
  }

  # 1. Indexed over cell type:
  for (u in 1:length(cts)) {
    suppressMessages(gc())

    # subset metadata and gene data by the barcodes for celltype u
    metsub = lme4metadata[lme4metadata$celltype == cts[u], ]
    df_ = as.matrix(gene_data[rownames(metsub), ])
    df_ = df_[match(rownames(df_), rownames(metsub)), ]

    if(is.null(celltype_genes_test)){
      gene_names = colnames(gene_data)
      print(paste0("fitting mixed model within ", print(cts[u])," for ",  length(gene_names), " genes" ))
    }
    if(!is.null(celltype_genes_test)){
      gene_names = celltype_genes_test[[u]]
      print(paste0("fitting mixed model within ", print(cts[u])," for ",  length(gene_names), " genes" ))
    }

    # 2. Indexed over gene
    for (i in 1:length(gene_names)) {
      # data for gene i
      dat_fit = cbind(metsub, gene = df_[ ,gene_names[i]])

      # fit model
      m1 = tryCatch(lme4::glmer(f1, data = dat_fit, family = poisson(link = "log")), error = function(e) return(NA))

      # return error cont. function
      if( suppressWarnings(is.na(m1))) {
        print(paste0(gene_names[i],  " could not fit  model for ", as.character(cts[u])))
        res_list[[i]] = NA
      } else{
        res = data.frame(t(summary(m1)$coefficients[2, ]))
        rnames = names(summary(m1)$coefficients[2, ])
        rnames[rnames == 'Estimate'] = 'fixed_effect'
        names(res) = rnames

        # return results
        res_list[[i]] = cbind(gene = gene_names[i],
                              res,
                              message = m1@optinfo$conv$lme4$messages,
                              model = f1)
      }
    }
    # combine fitted data for all genes
    res_list = res_list[!is.na(res_list)]
    resdf = do.call(rbind, res_list)
    res_celltype[[u]] = cbind(celltype = cts[u], resdf)
  }
  resdf_full = do.call(rbind, res_celltype)
  return(resdf_full)
}




