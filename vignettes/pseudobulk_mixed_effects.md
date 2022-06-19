Pseudobulk differential expression with nested group repeated measures
single cell experiment designs
================

-   <a href="#pseudobulk-mixed-effects-models"
    id="toc-pseudobulk-mixed-effects-models">1. pseudobulk mixed effects
    models</a>
-   <a href="#load-single-cell-data-aggregate-and-quaity-control"
    id="toc-load-single-cell-data-aggregate-and-quaity-control">Load single
    cell data, aggregate and quaity control</a>
-   <a
    href="#fit-models-and-specify-a-priori-contrasts-corresponding-to-the-desired-effect-comparisons"
    id="toc-fit-models-and-specify-a-priori-contrasts-corresponding-to-the-desired-effect-comparisons">Fit
    models and specify a priori contrasts corresponding to the desired
    effect comparisons</a>
-   <a
    href="#downstream-gene-set-enrichment-analysis-within-celltypes-for-different-effects"
    id="toc-downstream-gene-set-enrichment-analysis-within-celltypes-for-different-effects">Downstream
    gene set enrichment analysis within celltypes for different effects</a>
-   <a href="#further-analysis-and-curation-of-enrichment-results"
    id="toc-further-analysis-and-curation-of-enrichment-results">Further
    analysis and curation of enrichment results</a>
-   <a href="#extract-all-leading-edge-genes-indexed-by-cell-type"
    id="toc-extract-all-leading-edge-genes-indexed-by-cell-type">Extract all
    leading edge genes indexed by cell type</a>
-   <a
    href="#calculate-the-jaccard-similarity-of-the-leading-edge-genes-for-enrichments-within-a-given-cell-type-and-effect"
    id="toc-calculate-the-jaccard-similarity-of-the-leading-edge-genes-for-enrichments-within-a-given-cell-type-and-effect">Calculate
    the Jaccard similarity of the leading edge genes for enrichments within
    a given cell type and effect</a>
-   <a href="#visualization-of-enrichment-results"
    id="toc-visualization-of-enrichment-results">Visualization of enrichment
    results</a>

### 1. pseudobulk mixed effects models

These functions implement wrappers around limma for fitting fixed effect
linear models, and the `dream` method from the `variancePartition`
package forfitting mixed (i.e. varying) effects models. Mixed models are
necessary for experiments that have the design of repeated measurements
from the same donors in perturbation studies in multi sample ‘random’
(varying) effects. This is necessary to account for non-independence
when we have perturbation experiments with repeated measurements from
the same donors. To enable linear models (e.g. modeling the mean with a
normal distribution) to be fit to gene counts,
`variancePartition::dream` accounts for the mean variance trend via
incorporating voom observational weights. The fits are then shrunken
toward the genome wide trend using an empirical bayes step adjusting for
the per gene degrees of freedom estimated in the mixed model fits. This
approach is described in [Hoffman et al Bininformatics
2020](doi.org/10.1093/bioinformatics/btaa687) **this paper should be
cited if using the wrapper below `FitDream`.**

In this workflow, mixed effect model fits and gene set enrichment steps
are parallelized with the BiocParallel package.

``` r
#devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
suppressMessages(library(scglmmr))
suppressMessages(library(magrittr))
library(BiocParallel)
#the mixed model fits and gsea are parallelized 
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)
```

Below analysis of a 2 group repeated measures experiment design is
shown. Pre-treatment baseline effect differences, treatment effets
across all donors and the difference in treatment effects are all
compared while modeling variation in baseline expression using a random
intercept term.

| sample         |    sampleid    |           time | group          |      sex       |
|:---------------|:--------------:|---------------:|:---------------|:--------------:|
| 101_t0         |      101       |             d0 | low            |       F        |
| 101_t1         |      101       |             d1 | high           |       F        |
| 102_t0         |      102       |             d0 | low            |       M        |
| 102_t1         |      102       |             d1 | high           |       M        |
| … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) |

### Load single cell data, aggregate and quaity control

``` r

datapath = "mypath/"

# load seurat or sce object etc. 
s = readRDS("path/seuratobject.rds")

# define counts and metadata and subset to cells above rm seurat object from workspace 
meta = s@meta.data
umi = s@assays$RNA@counts
rm(s); gc()

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "celltype", sample_column = "sample")
tab$celltypes_remove; tab$`low representation celltypes`; tab$table

# remove cells prior to pseudobulk analysis 
meta = meta[!meta$celltype %in% tab$celltypes_remove, ]

# subset data 
umi = umi[ ,rownames(meta)]

# Create aggregated pseudobulk data
pb = scglmmr::PseudobulkList(
  rawcounts = umi,
  metadata = meta, 
  sample_col = "sample",
  celltype_col = "celltype",
  avg_or_sum = "sum"
  )


# Create aggretated sample level metadata
met = scglmmr::AggregateCellMetadata(
  cell.metadata = meta, 
  sample_column = 'sample', 
  variable_columns = c('subjectid', 'timepoint', 'response', 'age', 'sex'),
  pseudobulk.List = pb
  )

# creation of a combined grouping factor indicating both timepoint and group
# (e.g. t0_Group1, t0_group2) will make it simple to set up contrasts.
# here group means 'response group' i.e. high vs low responder. 
# also convert other variables to factors and do some standard transformation
# of metadata.
met$group.time = paste(met$group, met$timepoint, sep = '_')

# make sure this is a factor and the levels are in a coherent order for the 
# contrast matrix -- see below. here: 
# time 0 = 0, time 1 = 1. 
# low response = 0 high response = 1. 
met$group.time = factor(
  met$group.time,
  levels = c("1_0", "1_1", "0_0", "0_1")
  )

# now filter genes within each cell type that are reasonably expressed. 
design = model.matrix( ~ 0 + met$group.time)
dge = Normalize(pseudobulk.list = pb, design = design, minimum.gene.count = 5)
```

### Fit models and specify a priori contrasts corresponding to the desired effect comparisons

Below shows 2 steps. 1: Specify a varying intercept model using lme4
symbolic model formula. [More information on symbolic formula
representations of models](https://arxiv.org/pdf/1911.08628.pdfs). Most
experiment designs will have some covariate of interest such as time
relative to perturbation, other covariates which should be adjusted for
(such as batch or sex), and individuals measured with repeated
timepoints will have have a “subjectID” or “individual” variable
specifying which person the measurement came from. This variable will
often modeled using a varying intercept modeled with a normal
distribution which is specified with e.g. (1\| subjectID). Other formula
accepted by lme4 are also possible.

In a simple experiment, e.g. a one way anova design with baseline and
post treatment and no different outcome groups, if the time relative to
treatment is the target we want estimates for, we don’t need a design
matrix, we can simply extract the effect of “time” from the lme4 fits.

In the example below, we want estimate baseline differences, treatment
effects across all individuals and the difference in fold changes
between groups adjusting estimates for age and sex (and modeling the
individual variation with a varying effect). To do this we specify a
custom a. priori contrast matrix.

2: Using the factor variable combining group and time, we create
contrasts over the levels to define effects of interest. The function
`makeContrastsDream` is used for this purpose which works the same way
as the limma function makeContrasts. More simple contrasts or more
complex contrasts are also possible. This step is critical for deriving
the correct effec size estimates. More information on specifying
contrast matrices is available here: [A guide to creating design
matrices for gene expression
experiments](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)

``` r
# Now we  specify model 
f1 <- ~ 0 + age + sex + group.time + (1|subjectid) 

# make contrast matrix 
L2 = makeContrastsDream(
  formula = f1, 
  data = samplemd,
  contrasts = c(
    baseline = "group.time1_0 - group.time0_0", # note fixef model also fit for this single time point contrast 
    treatment_delta = "( group.time1_1 - group.time1_0 ) - ( group.time0_1 - group.time0_0 )",
    treatment = "( group.time1_1 + group.time0_1 ) / 2 - ( group.time1_0 + group.time0_0 ) / 2 "
    )
  )

# visualize the contrast matrix, compare to levels of group.time variable
# to check for correctly specified effects. 
plotContrasts(L2) 
```

Fit the models with `variancePartition::dream` which uses voom weights
in lme4 mixed effects models with an empirical bayes shrinkage toward
genome wide trend accounting for per gene degrees of freedom.

``` r

fit1 = FitDream(pb.list = dge, 
                sample.metadata = met, 
                lme4.formula = f1,
                dream.contrast.matrix = L2,
                ncores = 4)
```

Note, `variancePartition::dream` now incorporates an emperical Bayes
step derived for mixed effects models accounting for per gene degrees of
freedom, please see:
<https://github.com/GabrielHoffman/variancePartition/issues/54>.

### Downstream gene set enrichment analysis within celltypes for different effects

Run gene set enrichment analysis within each cell type fbased on genes
ranked by effect size for each of the effects defined above.

Here within each cell type and for each contrast we run gene set
enrichment analysis using 2 functions `ExtractResult` and `FgseaList`.
`ExtractResult` is used to extract a list of gene ranks (it can be used
as a wrapper around dream::topTable and limma::topTable to return a full
list of results). Since we had multiple covariates in a mixed model and
we used a custom contrast matrix, we specify the covariate of interest
from the contrast using arguments `coefficient.number` and `coef.name`
which are output in results based on the names of the contrasts.

With the list of ranks we then run gene set enrichment with the [fgsea
method by Korotkevich et
al](https://www.biorxiv.org/content/10.1101/060012v3) using the function
`FseaList`. Conveniently this also is parallelized using the same
biocparallel

Based on our contrast matrix (the object `L2` we created above with
`makeContrastsDream`), `coefficient.number = 1` corresponds to the
baseline difference between groups, `coefficient.number = 2` is the
difference in fold changes between the groups and
`coefficient.number = 3` is the pre vs post perturbation difference
across subjects in both groups. Estimates for the covariates are also
availabe.

#### Examples of extracting gene ranks based on contrasts and running gene set enrichment

``` r

# msigDB hallmark pathways are included in the scglmmr package
hlmk = scglmmr::hallmark
scglmmr::ExtractResult()
# extract genes ranked by treatment effect (across all donors) 
rtreat = ExtractResult(model.fit.list = fit1,
                       what = 'lmer.z.ranks',
                       coefficient.number = 3,
                       coef.name = 'treatment')
# run fgsea 
hlmk.treat = FgseaList(rank.list.celltype = rtreat,
                       pathways = hlmk,
                       BPPARAM = pparam)



# extract gene ranks of treatment effect (across all donors) 
rdelta = ExtractResult(model.fit.list = fit1,
                       what = 'lmer.z.ranks',
                       coefficient.number = 2,
                       coef.name = 'treatment_delta')
hlmk.delta = RunFgseaOnRankList(rank.list.celltype = rdelta,
                                pathways = hlmk,
                                BPPARAM = pparam)
```

On a technical note, it’s also possible to fit a simple model for the
baseline contrast not estimating the variation across subjects with a
random effect as this contrast is more typically estimated with a
standard group 1 vs 2 least squares approach. The function
`RunVoomLimma` is provided to fit these models. If you compare the
effect size estimates for an individual cell type for the
`coefficient.number = 1` from the mixed model fits (object `fit1` above)
to the models below they will be very highly corrleated along the
diagonal with some differences arising to the different degrees of
freedom, sample size and variance coming from the additional layer of
variation estimated by the mixed model. The choice of model to use is up
to the user.

``` r
# fit simple linear model for the baseline group level contrast 
design.2 = model.matrix(~0 + age + sex + group.time, data = met)
fit0 = scglmmr::RunVoomLimma(dgelists = dge, 
                           design_matrix = design.2, 
                           do_contrast_fit = T,
                           # we use only the first row of the contrast matrix L2
                           my_contrast_matrix = L2[ ,1])


# extract fixed efect model estimate of baseline contrast (pre treatment differences between groups)
r0 = ExtractResult(model.fit.list = fit0,
                   what = 'gene.t.ranks',
                   coefficient.number = 1,
                   coef.name = 'baseline')
## run gene set enrichment 
hlmk.0 = FgseaList(rank.list.celltype = r0,
                   pathways = hlmk,
                   BPPARAM = pparam)
```

### Further analysis and curation of enrichment results

``` r
# full set of leading edge genes indexed by celltype x effect x module 
lefull = scglmmr::GetLeadingEdgeFull(gsea.list = hlmk.treat,
                                     padj.filter = 0.02, 
                                     NES.filter = -Inf)

# extract model fit results instead of ranks
# automatically this sets argument `result` to 'statistics'
fit1.res = scglmmr::ExtractResult(model.fit.list = fit1 , 
                                  coefficient.number = 3, 
                                  coef.name = 'treatment')

# combine GSEA results with model coefficient for each gene in leading edge 
# include all leading edge genes irrespective of individual gene p value 
cr = scglmmr::CombineResults(gsealist = hlmk.treat, 
                             contrastlist = fit1.res, 
                             gseafdr = 0.02, 
                             genefdr = 1)
```

### Extract all leading edge genes indexed by cell type

This outputs a nested list of leading edge genes from each enrichment
across cell types. list level 1 is cell type, level 2 is enrichment
signal.

``` r
li = scglmmr::LeadingEdgeIndexed(gsea.result.list = hlmk.treat, padj.threshold = 0.02)
```

### Calculate the Jaccard similarity of the leading edge genes for enrichments within a given cell type and effect

This is useful for understanding which enrichment signals may come from
the same vs distinct genes.

``` r
# figpath.temp = here('figures')
treat.JI = EnrichmentJaccard(gsealist = hlmk.treat, 
                          indexedgenes = li, 
                          #saveplot = TRUE, 
                          #figpath = figpath.temp,
                          returnJaccardMtx = TRUE)

# curate results 
results.sorted = treat.JI$sortedgsea %>%
  dplyr::mutate(signal = paste(celltype, pathway, sep = '~'))
  
# data.table::fwrite(results.sorted, file = paste0(datapath, 'g0.result.sort.txt'),sep = "\t")
```

### Visualization of enrichment results

**Create a bubble plot heatmap of enrichment results within clusters**

``` r
# NES_filter to -Inf to plot positive and negative enrichment
p = PlotFgsea(gsea_result_list = hlmk.treat, NES_filter = -Inf,padj_filter = 0.05)
```

**Create heatmap of gene perturbation fold changes across cell subsets
based on model fit coefficients**

Extract and make a heatmap of log fold changes across celltypes.
HeatmapDiag uses the package
[slanter](https://CRAN.R-project.org/package=slanter) which accepts the
same arguments as `pheatmap`.

``` r
gene.mat = GetGeneMatrix(result.list = fit1, 
                         pvalfilter = 0.05, 
                         stat_for_matrix = 'logFC',
                         logfcfilter = 0.1)

#make heatmap  e.g. 
pheatmap::pheatmap(gene.mat, fontsize_row = 5)

# diagnoalize the heatmap so that the most correlated signals are grouped 
# this uses the package slanter 
HeatmapDiag(matrix = gene.mat, fontsize_row = 5)
```

**Create heatmap of gene perturbation fold changes from enrichments
across individuals within a cell subset**

Visualize the log counts per million at the sample level of the gene set
enrichment results.

``` r

# make tidy average data for visualization of weighted pb results 
lcmp = lapply(pb, edgeR::cpm, log = TRUE )

# 
le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = lcmp,
                                           gsea.list = hlmk.treat, 
                                           padj.filter = 0.1,
                                           NES.filter = -Inf)


# example plot of sample level average leading edge genes annotated 
scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr,
                               modulename = "HALLMARK_HYPOXIA",
                               elltype_plot = "CD4_NaiveTcell",
                               metadata = meta, 
                               metadata_annotate = c('group', 'timepoint', 'age', 'sex'),
                               sample_column = 'sample',
                               returnmat = FALSE, 
                               savepath = figpath, 
                               savename = "filename")

# see also: 
# scglmmr::TopGenesTidySampleExprs() - same as aobve for 'top' de genes estimated by model coeffcient / p value. 
# scglmmr::GetTidySummary() - for custom plotting


# sample level vsualization of average data above 
# repeat this for all enriched pathways 
heatpath = here("sime/path"); dir.create(heatpath)
for (i in 1:length(le_expr)) {
  cdat = le_expr[[i]]
  ctype = names(le_expr[i])
  umod = unique(cdat$module)
  for (u in 1:length(umod)) {
    scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr, 
                          modulename = umod[u], 
                          celltype_plot = ctype,
                          metadata = meta, metadata_annotate = c('group', 'timepoint', 'age', 'gender'),
                          sample_column = 'sample',
                          returnmat = F, 
                          savepath = heatpath,
                          savename = paste0(ctype, " ",umod[u],'.pdf'))
  }
}
```

You can also create a customized map by returning the matrix from
`LeadEdgeSampleHeatmap`.

``` r

# scglmmr function to extract leading edge genes 
lexp1 = LeadEdgeTidySampleExprs(av.exprs.list = lcpm, gsea.list = g1c, padj.filter = 0.05, NES.filter = 0)

# annotate time and batch on the heatmap
heatmap_anno = meta[, c('batch', 'timepoint')]
anno_color = list(
  timepoint = c('1' = "orange", '2' = 'red',  "0" = "white"),
  batch = c('1' = "black", '2' = "white")
)

# define your own custom color vector for the log fold change values
cu = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3", 
       "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")

# scglmmr function for leading edge gene matrix across donors
mat2 = scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr, 
                                 modulename = "reactome interferon signaling",
                                 celltype_plot = 'CD14_Mono',
                                 # metadata = meta, # unused  
                                 # metadata_annotate = c('batch'), # unused 
                                 sample_column = 'sample',
                                 returnmat = TRUE)
# draw heatmap 
pheatmap::pheatmap(mat2, 
                   border_color = NA,
                   treeheight_row = 0, treeheight_col = 10,
                   annotation = heatmap_anno,
                   annotation_colors = anno_color,
                   color = cu,
                   width = 5,  height = 7.6,
                   # this can be set to false to look at raw expression e.g. 
                   scale = "row",
                   filename = paste0(figpath, "mono_d1_ReactomeIFN.pdf")
                   )
```

``` r
sessionInfo()
```
