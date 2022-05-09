# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

#' SCMixedPoisson - Fit a gene level Poisson count mixed generalized linear model.
#'
#' @param gene_data matrix of raw single cell UMI counts genes as COLUMNS and cells as ROWS; i.e. a transpose t() of Bioconductor or Seurat raw.data slots if genes_test is specified, columns must contain all genes from genes_test i.e. umi_matrix = t(seurat_raw_counts[unique(unlist((genes_test)), ])
#' @param metadata metadata for each cell used for fitting model
#' @param model_formula example: 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)' the 'gene~' and '(1|subjectid)' variables should be kept the same. The variable 'timepoint' can be altered to be any treatment or perturbation effect being tested. User can also specify 'covariate_variables' to make formula gene~offset(nUMI) + covariate1 + covariate2 + test_variable + (1|subjectid) automatically.
#' @param test_variable the column of metadata coding the perturbation variable. e.g. 'timepoint' for formula: gene ~ timepoint + (1|subjectid)'
#' @param covariate_variables a vector of variables in metadata to add to the model as covariates; this is only used if model_formula is NULL; otherwise speify directly
#' @param celltype_genes_test An R list indexed by cell type: the subset of genes to test for each cell type.
#' @param save_path file path to save intermediate results for each cell type.
#'
#' @return a results matrix for all genes across all cell types tested.
#' @export
#'
#' @importFrom lme4 glmer
#' @importFrom emmeans emmeans
#' @importFrom stringr str_replace
#'
#' @examples
#' results = SCMixedPoisson(gene_data = gene_data,
#'   lme4metadata = meta_data,
#'   model_formula = 'gene ~ offset(log(nUMI)) + timepoint + (1|sampleid)',
#'   celltype_genes_test = celltype_indexed_gene_vector
#' )
SCMixedPoisson = function(gene_data,
                          metadata,
                          model_formula = 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)',
                          test_variable = NULL,
                          covariate_variables = NULL,
                          celltype_genes_test = NULL,
                          save_path){

  # specify and unify model formula and metadata for emmeans call
  if(!is.null(model_formula)){
    f1 = model_formula
  } else if(!is.null(covariate_variables)) {
    fixef = paste0(covariate_variables, collapse = " + ")
    f1 = paste("gene ~", fixef, '+', test_variable, '+ (1|subjectid)')
  }
  f1 = stringr::str_replace(string = f1, pattern = test_variable, replacement = 'perturbation')
  mform = as.formula(f1)
  names(metadata)[names(metadata) == test_variable] = "perturbation"
  ## Run pre model fit checks
  ctmd = split(metadata, f = metadata$celltype)
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
  if (is.null(test_variable)) {
    stop("specify metadata variable associated with pre-post perturbation, e.g. factor 'timepoint' w. levels c('0h', '24h')")
  }
  if (!c("celltype") %in% colnames(metadata)) {
    stop("in metadata the celltype / cluster column must be named 'celltype'")
  }
  if ( isFALSE(grepl("(1|subjectid)", f1, fixed = TRUE))) {
    stop("specify random effect (1|subjectid) in model_formula for the experimental units with repeated measurements (e.g. persons, mice)")
  }
  # Run checks for minimal required structure of f1 and f2 lme4 formulae
  if ( isFALSE(grepl("gene" , f1, fixed = TRUE))) {
    stop(" 'gene' must be specified in LHS of formula e.g. 'gene ~ ...'")
  }
  # init storage
  res_celltype = list()
  cts = as.character(unique(metadata$celltype))
  nct = length(cts)
  if(!is.null(celltype_genes_test)){
    stopifnot(all.equal(names(celltype_genes_test), cts))
  }

  # Indexed over celltype create data subset for celltype u
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    res_list = list() # clear the list for u + 1
    if(is.null(celltype_genes_test)){
      gene_names = as.character(unique(colnames(gene_data)))
      print(paste0("fitting mixed model within ", print(cts[u])," for ",  length(gene_names), " genes" ))
    }
    if(!is.null(celltype_genes_test)){
      gene_names = as.character(unique(celltype_genes_test[[u]]))
      ngene = length(gene_names)
      print(paste0("fitting models within ", print(cts[u])," for ", ngene, " genes" ))
    }

    # subset metadata and gene data by gene subset and cell barcodes for celltype u
    metsub = metadata[metadata$celltype == cts[u], ]
    df_ = as.matrix(gene_data[rownames(metsub), gene_names])
    df_ = df_[match(rownames(df_), rownames(metsub)), ]

    # 2) Indexed over genes within celltype u, run Poisson GLMM
    .dataenv <- environment()
    for (i in 1:length(gene_names)) {
      print(paste0('fitting ', gene_names[i], " ", i, ' of ', ngene,
                   ' genes for celltype ', u,' of ', nct, " ", cts[u]))

      dat_fit = cbind(metsub, gene = df_[ ,gene_names[i] ])
      dat_fit <- as.data.frame(dat_fit, env = .dataenv)

      m1 = tryCatch(
        lme4::glmer(formula = mform, data = dat_fit, family = poisson(link = "log")),
        error = function(e) return(NA)
      )
      emm1 = tryCatch(
        emmeans::emmeans(object = m1, specs = revpairwise ~ perturbation)$contrast,
        error = function(e) return(NA)
      )
      # error handling on model fit and marginal means
      if(suppressWarnings(is.na(m1))) {
        res_list[[i]] = NA
      }
      if(suppressWarnings(is.na(emm1))) {
        res_list[[i]] = NA
      }
      else{
        # extract convergence error messsages and enforce same data frame structure on messages
        messages = m1@optinfo$conv$lme4$messages
        mdf = data.frame()
        if (is.null(messages)) {
          mdf = data.frame('message 1' = 'none', 'message 2' = 'none')
        }
        # note genes with singular fits
        if (length(messages) == 1){
          mdf = data.frame('message 1'= as.character(m1@optinfo$conv$lme4$messages[[1]]),
                           'message 2' = 'none',
                           stringsAsFactors = FALSE)
        }
        # Flag genes with convergence errors related to Hessian gradient
        if (length(messages) > 1) {
          mdf = data.frame(
            'message 1'= as.character(m1@optinfo$conv$lme4$messages[[1]]),
            'message 2' = as.character(m1@optinfo$conv$lme4$messages[[2]]),
            stringsAsFactors = FALSE
          )
        }
        # extract results onto a single line
        res_list[[i]] = cbind.data.frame(
          'gene' = as.character(gene_names[i]),
          'celltype' = as.character(cts[u]),
          as.data.frame(emm1[1, ]),
          'model' = f1,
          mdf,
          stringsAsFactors = FALSE
        )
      }
    }
    # rbind fitted model results for each gene within celltype u
    res_list = res_list[!is.na(res_list)]
    res_list = res_list[!is.null(res_list)]
    d = as.data.frame(do.call(rbind, res_list))
    # calculate adjusted p values
    d$contrast_padj = p.adjust(p = d$p.value, method = "BH")
    saveRDS(object = d, file = paste0(save_path, cts[u], 'result.rds'))
    # store for celltype u
    res_celltype[[u]] = d
  }
  # rbind fitted model results for all u celltypes
  resdf_full = as.data.frame(do.call(rbind, res_celltype))
  saveRDS(object = resdf_full, file = paste0(save_path, 'SCMixedPoisson_full_result.rds'))
  return(resdf_full)
}


