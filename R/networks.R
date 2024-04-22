

#' Jmatlist
#'
#' @param indexedgenes returned by scglmmr::LeadingEdgeIndexed
#' @return a list indexed by cell type of the jaccard index within each cell type of the leading edge gene content from each enrichment.
#' @importFrom GeneOverlap newGOM getMatrix
#' @export
#'
#' @examples
#'\dontrun{
#'# EXAMPLE:
#'gsea.list = FgseaList()
# indexed.gene.list = LeadingEdgeIndexed(gsea.result.list = gsea.list, padj.threshold = 0.05)
# jmatrix = Jmatlist(indexedgenes = indexed.gene.list)
#' }
Jmatlist = function (indexedgenes, ...) {
  jmat = list()
  nenriched = unlist(lapply(indexedgenes, length))
  # store NA for list length 0 if user does not filter list first.
  for (i in 1:length(indexedgenes)) {
    print(i)
    if (nenriched[i] == 0) {
      jmat[[i]] = NULL
    }
    if (nenriched[i] == 1) {
      mod = names(indexedgenes[[i]])
      diagmat = matrix(1, dimnames =  list(
        paste(names(indexedgenes)[i], mod, sep = ': '),
        paste(names(indexedgenes)[i], mod, sep = ': '))
      )
      jmat[[i]] = diagmat
    }
    if (nenriched[i] > 1) {
      # compute jaccard index for ech module pair.
      overlap = GeneOverlap::newGOM(gsetA = indexedgenes[[i]],
                                    gsetB = indexedgenes[[i]],
                                    genome.size = NULL)
      jaccard_matrix = GeneOverlap::getMatrix(object = overlap, name = "Jaccard")
      colnames(jaccard_matrix) = paste(names(indexedgenes)[i], colnames(jaccard_matrix), sep = ': ')
      rownames(jaccard_matrix) = paste(names(indexedgenes)[i], rownames(jaccard_matrix), sep = ': ')
      jmat[[i]] = jaccard_matrix
    }
  }
  names(jmat) = names(indexedgenes)
  # remove cell types with no enrichments
  jmat = base::Filter(x = jmat, f = length)
  return(jmat)
}



#' Signaturescores
#'
#' @param feature.list pseudobulk list indexec by cell type (the list of DGEList objects)
#' @param indexedgenes object returned by scglmmr::LeadingEdgeIndexed
#' @return a list indexed by cell type of  sample level scores (average module Z score as in Kotliarov et al Nature Medicine 2020) of each signature leading edge genes derived from statistical models.
#' @importFrom data.table rbindlist
#' @export
#'
#' @examples
#'\dontrun{
#'
#'index = LeadingEdgeIndexed(gsea.result.list = glist,padj.threshold = 0.05)
#'jmat = Jmatlist(index)
#'sig = Signaturescores(feature.list = gexp, indexedgenes = index)
#'
#' }
#'
#'

Signaturescores = function(feature.list, indexedgenes, ...) {

  print("checking list indices are identical");print(names(indexedgenes)); print(names(feature.list))

  # fix this for the user if not ordered properly
  if (!isTRUE(all.equal( names(indexedgenes), names(feature.list)))) {
    warning(' : filtering and reordering `gene.expr` by the indices of `indexedgenes` test  names(indexedgenes); names(gene.expr)')
    feature.list = feature.list[names(indexedgenes)]
    print('reordered input data: ');print(names(indexedgenes));print(names(feature.list))
  }
  # check
  stopifnot(isTRUE(all.equal( names(indexedgenes), names(feature.list))))

  nenriched = unlist(lapply(indexedgenes, length))
  av.zscore = list()
  for(i in 1:length(feature.list)) {
    if(isTRUE(nenriched[i] == 0)){
      av.zscore[[i]] = NA
    }
    if(isTRUE(nenriched[i] > 0)){
      av = base::as.matrix(feature.list[[i]])
      res = list()
      for (u in 1:length(indexedgenes[[i]])) {
        mod.genes = indexedgenes[[i]][[u]]

        # test signal for non overlapping genes
        test.genes = length(mod.genes) - length(mod.genes)
        if (test.genes > 0) {
          warning('check module genes in indexedgenes object - repeated gene symbols within the same signal detected ')
        }

        mod.genes = intersect(unique(mod.genes), rownames(av)) # note mod.genes should already be unique.
        avz = t(scale(t(av[mod.genes, ])))
        xz = as.data.frame(t(colMeans(avz, na.rm = T)))
        res[[u]] = xz
      }
      res = data.table::rbindlist(l = res,fill = TRUE)
      res = as.data.frame(res)
      rownames(res) = paste(names(feature.list)[i], names(indexedgenes[[i]]),sep = ': ')
      # store
      av.zscore[[i]] = res
    }
  }
  # name the results by cell type
  names(av.zscore) = names(feature.list)
  return(av.zscore)
}








#' SLImatrix
#'
#' @param sig.scores returned by scglmmr::EnrichmentJaccard
#' @param jmat.list, returned by scglmmr::Jmatlist
#' @param correlation.type string specifying correlation type as calculated by Hmisc::rcorr either default 'spearman' or 'pearson'
#' @return a list indexed by cell type. first element of the list is the sli.adjusted corelation matrix (within cell type corerlations are adjusted for shared gene content -- see Mulè et al 2024 Immunity (in press) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10120791/), second is the original object returned by Hmisc::rcorr.
#' @importFrom Hmisc rcorr
#' @export
#'
#' @examples
#'\dontrun{
#'index = LeadingEdgeIndexed(gsea.result.list = glist,padj.threshold = 0.05)
#'jmat = Jmatlist(index)
#'sig = Signaturescores(feature.list = gexp, indexedgenes = index)
#'sli.list = SLImatrix(sig.scores = sig, jmat.list = jmat, correlation.type = 'spearman')
#'
#' }
#'

SLImatrix = function(sig.scores, jmat.list, correlation.type = 'spearman'){

  stopifnot(isTRUE(all.equal(names(jmat.list), names(sig.scores))))
  # combine and calculate correlation across full matrix
  dm = bind_rows(sig.scores)
  cor.test = Hmisc::rcorr(t(as.matrix(dm)),type = correlation.type)
  # unadjusted rho
  # intracellular correlations
  ds.cor = lapply(sig.scores, function(x) Hmisc::rcorr(t(x), type = correlation.type)$r)

  sli = list()
  for (i in 1:length(jmat.list)) {
    stopifnot(isTRUE(all.equal(rownames(ds.cor[[i]]), rownames(jmat.list[[i]]))) )
    sli[[i]] = ds.cor[[i]] - jmat.list[[i]]
  }

  # replace the matrix values from intracellular correlations with the SLI values
  mat = cor.test$r
  for (i in 1:length(sli)) {
    print(i)
    # get index of cols and rows
    row.replace = which(rownames(mat) %in% rownames(sli[[i]]))
    col.replace = which(colnames(mat) %in% colnames(sli[[i]]))

    # replace values along the square diagonal of the matrix
    stopifnot(isTRUE(all.equal(row.replace, col.replace)))

    # check structure
    # stopifnot(isTRUE(all.equal(mat[row.replace, col.replace],ds.cor[[i]])))
    # full spearman matrix subset by rows of celltype i
    # original spearman correlation matrix for celltype i
    # replace iteratively
    mat[row.replace,col.replace] = sli[[i]]
    # this should be going down with each iteration
    print(sum(mat))
  }
  # celltypes with only 1 enrichment intracellular not replaced iteratively; replace here.
  diag(mat) = 0
  # print(mat %>% pheatmap::pheatmap(cluster_cols = F, cluster_rows = F, main = 'corrected'))
  ret.list = list('sli.adjusted.matrix' = mat, 'original.correlation.object' = cor.test)
  return(ret.list)
}



#' p.adjust.cormat
#'
#' @param hmisc.cor runs on individual objects returned by hmisc. returned by scglmmr::EnrichmentJaccard
#' @param method argument to stats::p.adjust - defaults to 'fdr'
#' @return adjusted p value matrix corresponding to correlation matrix
#' @importFrom Matrix::isSymmetric
#' @export
#'
#' @examples
#'\dontrun{
#'index = LeadingEdgeIndexed(gsea.result.list = glist,padj.threshold = 0.05)
#'jmat = Jmatlist(index)
#'sig = Signaturescores(feature.list = gexp, indexedgenes = index)
#'sli.list = SLImatrix(sig.scores = sig, jmat.list = jmat, correlation.type = 'spearman')
#'
#'fdr.p = list()
#'for (i in 1:length(sli.list)) {
#'  mat2 = sli.list[[i]]$sli.adjusted.matrix
#'  fdr.p[[i]] = p.adjust.cormat(hmisc.cor = sli.list[[i]]$original.correlation.object, method = 'fdr')
#' }
#'
p.adjust.cormat = function(hmisc.cor, method = 'fdr'){
  stopifnot(isTRUE(Matrix::isSymmetric(hmisc.cor$P)))
  p.adj.lower = stats::p.adjust(hmisc.cor$P[lower.tri(hmisc.cor$P)], method = method)
  p.adj.upper = stats::p.adjust(hmisc.cor$P[upper.tri(hmisc.cor$P)], method = method)
  p.adj.mx <- matrix(rep(0,ncol(hmisc.cor$P)*ncol(hmisc.cor$P)), nrow = ncol(hmisc.cor$P))
  p.adj.mx[lower.tri(p.adj.mx)] <- p.adj.lower
  p.adj.mx[upper.tri(p.adj.mx)] <- p.adj.upper
  diag(p.adj.mx) = 1
  colnames(p.adj.mx) = rownames(p.adj.mx) = colnames(hmisc.cor$P)
  stopifnot(isTRUE(Matrix::isSymmetric(p.adj.mx)))
  return(p.adj.mx)
}



