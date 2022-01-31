scglmmr *S*ample-level *S*ingle-cell *G*eneralized *L*inear *M*ultilevel
*M*odels in *R*
================
Matt Mulè

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

``` r
devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
library(scglmmr)
```

## readme updates TO DO:

add TOC  
add makeContrastsDream  
change the workflow to represent the most up to dat PB workflow  
Add Enrichmentjaccard  
add new single cell models module and gene
level

## Single cell within cluster perturbation response differential expression

The purpose of this software is to analyze single cell genomics data
with pre and post perturbation measurements from the same individuals,
including complex designs wehre individuals with repeated measurements
are nested within in multiple response groups. The focus is on
implementing flexible generalized linear multilevel models to derive
group (i.e. good or poor clinical outcome, high or low rug response) and
treatment associated effects *within cell types* defined either by
protein (e.g. with CITE-seq data) or transcriptome based clustering
followed by downstream enrichment testing and visualization.

By default, the effect of treatment/perturbation across all subjects,
the baseline differences between outcome groups, and the difference in
the treatment effect between the outcome groups are tested. Any number
of model covariates can be specified and by default the package uses a
random intercept model to accomodate the non-independence of expression
within each subject.

An overview of scglmmr methods:

### 1\. pseudobulk mixed effects models

These functions implement convenience wrappers around both `limma` for
fitting traditional linear models (for e.g. baseline differenced between
groups) and the `dream` method from the `variancePartition` package.
This is the only pseudobulk differential expression method that can
accomodate both fixed and random effects (which are statistically
necessary to account for non-independence of multiple measurements from
the same donor) via use of `lme4` for mixed effects modeling while
accounting for measurement uncertainty via pooling unqeual library sizes
by incorporating `voom` observational weights. This method performed
well in related simulation studies by [Crowell et.
al. 2020](https://www.nature.com/articles/s41467-020-19894-4). scglmmr
performs statistical contrasts between groups comparing the least
squares means after accounting for model covariates is empowered by the
[emmeans](https://cran.r-project.org/web/packages/emmeans/index.html)
package.

### 2\. single cell gene level mixed effects models

A function for implemeting single cell gene-level poisson GLMM within
clusters.

### 3\. single cell module level mixed effect models

Functions for calculating single cell weighted and non weighted module
scores with mixed effects models fit on these module scores.

### 4\. Downstream enrichment testing and visualization

scglmmr implements enrichment testing of pseudobulk or single cell
results with wrappers around methods from the [fast set gene enrichment
(fgsea)](https://www.biorxiv.org/content/10.1101/060012v2#disqus_thread)
and
[clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
R packages.

### Philosophy

The scglmmr package considers each cluster/ cell type as a separate
‘experiment’ similar to fitting separate models to different FACS
sorted leukocyte subsets followed by RNAseq. Clusters can be defined by
transcriptome clustering and the approach is particularly well suited
for CITE-seq data with cells clustered based on normalized protein
expression levels which can be improved by normalizing and denoising
with our method, [dsb](https://github.com/niaid/dsb).

**Experiment designs (within each cluster / celltype) supported by
scglmmr** Below is a 2 group repeated measures experiment. A random
intercept model is required to account for autocorrelation in each
donors post perturbation response with their baseline level of
expression. This data can be accomodated by scglmmr. More simple
experiment designs are also supported, for example data with 2 groups
but not repeated pre/post treatment measurements, see the
`RunVoomLimma()`
function.

| sample         |    sampleid    |      timepoint | Group          |      sex       |
| :------------- | :------------: | -------------: | :------------- | :------------: |
| 101\_t0        |      101       |             d0 | good           |       F        |
| 101\_t1        |      101       |             d1 | poor           |       F        |
| 102\_t0        |      102       |             d0 | good           |       M        |
| 102\_t1        |      102       |             d1 | poor           |       M        |
| … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) |

### 1\. “Pseudobulk” mixed effects model

In the case of pseudobulk libraries derived from single cell RNAseq
data, interpolated model weights should correlate with library sizes
within each cell type. Indeed we confirmed that the log count per
million and square root standard deviation of genes had the expected
monotonically decreasing trend within each cell type and that the
resulting interpolated model weights across genes in a given sample were
highly correlated with both the number of cells used to create the
library and the sample’s total mRNA library size. Later, a more recent
report introduced the muscat R package and simulation studies for
testing within cluster differential expression. 10.1101/713412. This
paper included simulation and validation studies of established bulk
RNAseq methods for within cluster differential expression analysis. In
that report, a pseudobulk method that performed particularly well was
using limma with voom observational weights, and mixed effects (i.e. the
ability to specify a random effect for donor baseline expression) with
lme4 model formula implemented through the “dream” method available
through the vairance partition package. We implemented a weighted mixed
effects model using this limma + voom + lme4 method within each cell
type using a random intercept for each donor with contrast coding to
test for group level fold change differences.

``` r
#devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
library(scglmmr)
datapath = "mypath/"

# load seurat or sce object etc. 
s = readRDS("path/seuratobject.rds")

# define counts and metadata and subset to cells above rm seurat object from workspace 
meta = s@meta.data
umi = s@assays$RNA@counts
rm(s); gc()

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "lineage", sample_column = "sample")
tab$celltypes_remove; tab$`low representation celltypes`; tab$table

# remove cells prior to pseudobulk analysis 
meta = meta[!meta$celltype_label_3 %in% tab$celltypes_remove, ]

# subset data 
umi = umi[ ,rownames(meta)]

# pseudobulk workflow 
pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample",
                             celltype_col = "lineage", avg_or_sum = "sum")
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",
                                      variable_column = "cohort_timepoint", pseudobulklist = pb)
dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)

# custom a priori contrasts over a combined factor of the group and timepoint variables 
c_mat = makeContrasts(
  foldchange_difference = (group_timepoint1_1 - group_timepoint1_0) - (group_timepoint0_1 - group_timepoint0_0),
  time1_foldchange = (group_timepoint1_1 + group_timepoint0_1) / 2  - (group_timepoint1_0 + group_timepoint0_0) / 2,
  baseline_groups = (group_timepoint1_0 - group_timepoint0_0),
  levels = colnames(designmat)
)

#### REPLACE THE ABOVE WITH makeContrastsDream see v3aso3 contrast model in fsc. 


  
# fit mixed model for the multi timepoint contrasts 
fit = scglmmr::dreamMixedModel(dge_lists = dge, 
                               apriori_contrasts = TRUE, 
                               sample_column = 'sample',
                               cell_metadata = meta, 
                               contrast_matrix = c_mat, 
                               design_matrix = designmat, 
                               lme4_formula =  '~ 0 + age + gender + cohort_timepoint + (1|sampleid)', 
                               fixed_effects = c('age', 'gender', 'cohort_timepoint'), 
                               plotsavepath = figpath, 
                               ncores = 4)

# fit simple linear model for the baseline group level contrast 
bl = scglmmr::RunVoomLimma(dgelists = dge, 
                           design_matrix = designmat, 
                           do_contrast_fit = T,
                           my_contrast_matrix = c_mat[ ,3])

# save 
saveRDS(bl, file = paste0(datapath, "blfit_object.rds"))
saveRDS(fit, file = paste0(datapath, "fit_object.rds"))
```

### Downstream analysis on pseudobulk results

Here we provide examples of downstream analysis on results of the model
fit above focusing on the first contrast (the first coefficient in the
results), corresponding to:  
foldchange\_difference = (group\_timepoint1\_1 - group\_timepoint1\_0) -
(group\_timepoint0\_1 - group\_timepoint0\_0)

``` r
# hlmk = readRDS(file = here("signature_curation/hallmark.rds"))
figpath = "your/path"
test = scglmmr::GetRankResults(limma.fit.object.list = bl, coefficient.number = 1, "test")
res = scglmmr::GetContrastResults(limma.fit.object.list = bl, coefficient.number = 1, contrast.name = "test")

# get the mixed effects model results for the difference in fold changes 
fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit, 
                                         coefficient.number = 1,
                                         contrast.name = "foldchangedifference")
fit_rank = scglmmr::GetRankResultsRaw(contrast.result.raw.list = fit_res)


# Gene Set Enrichment Analysis 
gsea1 = scglmmr::RunFgseaOnRankList(rank.list.celltype = test, pathways = hlmk)
d = scglmmr::RbindGseaResultList(gsea_result_list = gsea1,NES_filter = -Inf,padj_filter = 0.2)
scglmmr::GSEABubblePlot(d, save_path = figpath, save_name = "plot.pdf")
# also see scglmmr::GSEABarPlot()

# full set of leading edge genes vs celltypes 
lefull = scglmmr::GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.1,NES.filter = -Inf)

# combine GSEA results with model coefficient for each gene in leading edge 
cr = scglmmr::CombineResults(gsealist = gsea1, contrastlist = res, gseafdr = 0.1, genefdr = 1)


# make tidy average data for visualization of weighted pb results 
av = scglmmr::PseudobulkList(rawcounts = umi, 
                             metadata = meta, 
                             sample_col = "sample", 
                             celltype_col = "celltype",
                             avg_or_sum = 'average')
le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av,
                                           gsea.list = hlmk_ctm0, 
                                           padj.filter = 0.1,
                                           NES.filter = -Inf)


# example plot of sample level average leading edge genes annotated 
scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr,
                               modulename = "MODULENAME",
                               elltype_plot = "TCELL",
                               metadata = meta, 
                               metadata_annotate = c('group', 'timepoint', 'age', 'gender'),
                               sample_column = 'sample',
                               returnmat = F, 
                               savepath = figpath, 
                               savename = "filename")

# see also: 
# scglmmr::TopGenesTidySampleExprs() - same as aobve for 'top' de genes estimated by model coeffcient / p value. 
# scglmmr::GetTidySummary() - for custom plotting


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

# plot only the leading edge genes from enrichment in a heatmap of log Fold Change estimated for a model contrast 
le = scglmmr::GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.1, NES.filter = -Inf)
genesub = do.call(rbind, le) %$% gene %>% unique 
mtx2 = scglmmr::GetGeneMatrix(result.list = res, 
                    stat_for_matrix = "logFC",
                    gene_subset = genesub, 
                    pvalfilter = -Inf, 
                    logfcfilter = 0.1)

pheatmap::pheatmap(mtx2,
                   breaks =seq(from = 0, to = 2,length.out = 99),
                   filename = paste0(figpath,"LEgenes_heatmap.pdf"))



### Hypergeometric erichment 
load(termdf) # this term2gene dataframe is included in the package see clusterProfiler
hyp = scglmmr::RunHypergeometricTest(result_list = fit_res,
                                     TERM2GENE_dataframe = termdf,
                                     pval_threshold = 0.1,
                                     logFC_threshold = 0,
                                     usefdr_threshold = FALSE)
# plot results 
scglmmr::PlotHypergeometric(hyperg_result = hyp, 
                            p.adjust.filter = 0.1,
                            genenumber_filter = 2,
                            savepath = figpath,
                            savename = "name", 
                            title = "title")

# calculate average module z score across samples (z score across samples of gene z scores)
# baseline_samples = c('') # user defined vector 
lapply(av, function(x), x[ ,baseline_samples])
mz = scglmmr::AverageSampleModuleZscore(average.metacell.list = av,
                                        module.list = hlmk,use.module.subset = F)
```

### 2\. Single cell gene level DE testing

This single cell level models of differential expression pre and post
treatment testing of genes by fitting a Poisson mixed glm. You specify a
reduced and full model formula to extract p values and effect sizes.

Recommended: The function accepts an argument ‘celltype\_genes\_test’
this is a named list indexed by celltypes.

First we load a single cell object an define highly variable genes
within each celltype, this can be done many different ways.

``` r
# Single cell mixed model. 
library(scglmmr)
suppressMessages(library(magrittr))
set.seed(1990)

datapath = 'your_path_to_save_results'
s = readRDS('path_to_seurat_or_SingleCellExperimentObject')

# define highly variable genes within each cell type 
gene_list = list()
for (i in 1:length(unique(s@meta.data$celltype))) {
  gene_list[[i]] = subset(s, celltype_joint == celltype_unique[i]) %>% 
    Seurat::FindVariableFeatures(x[[i]],top.genes = 10000)
  gene_list[[i]] = Seurat::VariableGenes(gene_list[[i]])
}

# one can also use scran trendVar etc. 
```

Now we run the single cell Poisson model by specifying a model formula,
the formula must have the LHS equal to `gene ~` and include a random
effect term `(1|subjectid)`. the variable that corresponds to the
perturbation variable as a 2 level factor with pre and post perturbation
measurements (repeated measurements for the same subject denoted by
`subjectid` must be identified, for example, if the metadata variable
coding for the time post perturbation is called treatment.time, then the
model formula could be `f1 = 'gene ~ offset(log(nUMI)) + treatment.time
+ (1|subjectid)'` with this model one would set the argument
`test_variable = treatment.time` which should have 2 ordered levels
corresponding to pre and post treatment. Within the function, the
`emmeans` package automatically tests the difference in the marginal
means of the post vs pre perturbation effect.

``` r
# load subset of genes to test for each celltype 
gene_union = unique(unlist(genes_test))

# create day 1 metadata dataframe from a single cell object (seurat, SingleCellExperiment etc.)
s@meta.data$subjectid = factor(as.character(s@meta.data$sampleid))

# create metadata for lme4
md = s@meta.data %>%
  droplevels() %>%
  dplyr::select(barcode_check, batch, subjectid, timepoint, celltype_joint, nUMI) %>%
  mutate(timepoint = factor(timepoint, levels = c("d0", "d1"))) %>% 
  rename(celltype = celltype_joint) %>% 
  column_to_rownames('barcode_check')

# get gene data 
gene_data = Matrix::t(s@raw.data[gene_union, ])

# we can now remove whatever single cell object was loaded in the environment. 
rm(s); gc()

## specify model formula 
f1 = 'gene ~ offset(log(nUMI)) + timepoint + (1|subjectid)'

# method 1
results = scglmmr::SCMixedPoisson(gene_data = gene_data,
                                  metadata = md,
                                  model_formula = f1, 
                                  test_variable = 'timepoint',
                                  celltype_genes_test = genes_test,
                                  save_path = datapath)
# write results 
# data.table::fwrite(x = results,file = paste0(datapath, "/h1_sc_time_merged_v2.txt"), sep = '\t')
```

Another way to specify the model automatically:

``` r
results = scglmmr::SCMixedPoisson(gene_data = gene_data,
                                  metadata = md,
                                  test_variable = 'timepoint',
                                  covariate_variables = c('BMI', 'SEX'),
                                  celltype_genes_test = genes_test, 
                                  save_path = datapath)
```

This would automatically test the delta for pre and post perturbation
calculating the marginal means over the coavariates BMI and SEX with the
generalized linear mixed effects Poisson model: `'gene ~ BMI + SEX +
timepoint +
(1|subjectid)'`.

### enrichment of single cell poisson model of timepoint effect using clusterprofiler

``` r
## enrichment on these results 

### plot enrichment with a hypergeometric test against the li BTM 
scd1 =  split( results , f = mmres$celltype )
scd1 = lapply(scd1, function(x){x %>%
    dplyr::select(gene, celltype, estimate, padj)  %>%
    dplyr::rename('logFC'= fixed_effect,  'P.Value' = padj) })


load(termdf)
hypsc = scglmmr::RunHypergeometricTest(result_list = scd1, 
                                       TERM2GENE_dataframe = term_df,
                                       pval_threshold = 0.05,
                                       logFC_threshold = 0.1)  
scglmmr::PlotHypergeometric(hyperg_result = hypsc,
                            p.adjust.filter = 0.05,
                            genenumber_filter = 2,
                            savepath = figpath, 
                            savename = "d1singlecell_hypergeometric_padj0.05.pdf", 
                            title = "single cell 24h vs baseline gene_exprs ~ offset(nUMI) + timepoint + (1|subjectid) \n 
                            hypergeometric test LI BTM, FDR 0.05, gene filter n > 2", 
                            height = 5.5, width = 8)
```

### 2\. Single cell module level testing + group level fold change comparisons

This is designed for testing the difference in single cell module
activity scores between groups at baseline, between groups, post
treatment across groups and the difference in treatment effect between
groups using the same *a priori* contrast approach as in the pseudobulk
method described above.

#### Part I

You first fit a module score to each single cell. Options: *threshold*
what percent of genes in the module must a cell express to get a score
(i.e. of 0.1, if the cell has less than 10% of genes \>0 in the module
don’t score) *return\_weighted* whether to return the weighted module
score which is the average \* weight where weight = number of genes in
the cell in that module with non-zeroexpression

*cellwise\_scaling* should the scores be scaled across cells? Not
recommended unless you are testing only a single celltype.
*Seurat\_version* if this is set to “2” it will accomodate seurat v2, if
set to “3” will accomodate v3

#### Part II fitting models:

*The SCGroupContrastGLMM Function* To simultaneously test the difference
in the fold changes due to treatment (for example) the levels of the
combined group + time factor "group\_id have to be ordered as follows
(can be any names, these are just the order): group1 timepoint 0, group1
timepoint 1, group2 timepoint 0, group2 timepoint 1

It is also possible to use custom contrasts by changing the argument to
contrast\_list from contrast\_2 (automatically loaded) to something
else.

Specify any lme4 model formula in the argument to f1, with the
randomeffect term for subjectid set to `(1|sampleid)` which should have
a pre and post perturbation measurement. `f1 = 'modulescore ~ age +
gender + group_id + (1|sampleid)'`

This function will also automatically save a helpful plot of the least
squares means (the residual means of each group after accounting for
covariates in the model) of the module over the levels of each combined
group\_timepoint factor.

``` r

# load data
Seurat = readRDS("my_seurat_object.rds")

# add cellwise module score for each signature 
mod_scores = WeightedCellModuleScore(seurat_object = Seurat,
                                     module_list = btm,
                                     threshold = 0.1,
                                     return_weighted = FALSE, cellwise_scaling = FALSE, 
                                     Seurat_version = "2") 
Seurat = AddMetaData(Seurat,metadata = mod_scores)
module_n = names(sig_test)

# set up module data frame 
module_df = Seurat@meta.data %>% select(barcode_check, celltype_joint, module_n) 

# format metadata as factors group_id is order leveled for contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
md = Seurat@meta.data %>% 
  mutate(group_id = factor(treat_time,  levels = c('pre_low', 'post_low', 'pre_high', 'post_high'))) %>%
  mutate(sampleid = factor(sampleid)) %>% 
  select(barcode_check, celltype_joint, sampleid,  age, group_id)

# Fit mixed model 
plot_savepath = paste0(my_figure_save_path, "/marginalmeans/"); dir.create(plot_savepath)

# specify any random intercept model e.g.
f1 = 'modulescore ~ age + group_id + (1|sampleid)'

# fit sc mod mixed model on ewighted module scores. 
mm_res = SCGroupContrastGLMM(module_data_frame = module_df, 
                             celltype_column = 'celltype',
                             metadata = md,
                             fixed_effects = NULL,
                             lmer_formula = f1,
                             plotdatqc = TRUE,
                             figpath = 'your/file/path')

### This function returns p values, effect sizes for baseline and difference in fld changes between the groups 
colnames(mmres)
# "celltype"                    "modulename"                  "contrasttime1vs0_group2vs1" 
# "estimatetime1vs0_group2vs1"  "std.errortime1vs0_group2vs1" "dftime1vs0_group2vs1"       
# "z.ratiotime1vs0_group2vs1"   "p.valuetime1vs0_group2vs1"   "contrasttime0_group2vs1"    
# "estimatetime0_group2vs1"     "std.errortime0_group2vs1"    "dftime0_group2vs1"          
# "z.ratiotime0_group2vs1"      "p.valuetime0_group2vs1"      "d0 low_marginal_mean"       
# "d1 low_marginal_mean"        "d0 high_marginal_mean"       "d1 high_marginal_mean"      
# "formula"                     "statistictime1vs0_group2vs1" "statistictime0_group2vs1"   
```

A plot like this will be generated for each module for each cell type
with the estimated marginal means in the right margin and the actual
single cell module score distribution across each level of
group\_id:

![image](https://user-images.githubusercontent.com/15280712/89591805-1f5a1200-d819-11ea-8b49-5b84477dd178.png)

*The WilcoxWithinCluster function* To run a simple test of baseline
differences between groups, use the *WilcoxWithinCluster* function as
below if you only have one timepoint but 2 outcome groups.

``` r

######## baseline differences with a simple wilcox test 
weighted_score = WeightedCellModuleScore(seurat_object = s, module_list = sig,
                                         threshold = 0.1, cellwise_scaling = FALSE,
                                         return_weighted = TRUE,
                                         Seurat_version = "3")
 
# set module dataframe
module_df = cbind(celltype_joint = as.character(s@meta.data$seurat_clusters),
                  barcode_check = s@meta.data$barcode_full,
                  weighted_score) %>% 
  droplevels()
 
# format lme4 metadata (note here just wilcox so group IDS do not have to be lme4 factors
meta_data = s@meta.data %>%
  droplevels() %>%
  select(celltype_joint = seurat_clusters,  barcode_check = barcode_full, sample, group_id = IRAE) %>%
  filter(group_id %in% c("poor_outcome", "good_outcome")) %>%
  mutate(group_id = factor(group_id, levels = c("poor_outcome", "good_outcome")))
 
# Fit wilcox test and do bonferonni correction 
plot_savepath2 = file.path(figpath, "module_distribution/"); dir.create(plot_savepath2)
m2 = module_df %>% filter(barcode_check %in% meta_data$barcode_check )
cvnc = WilcoxWithinCluster(module_data_frame = m2, lme4metadata = meta_data, plot_savepath = plot_savepath2)
cvnc$padj = p.adjust(p = cvnc$p.value.t0.wilcox, method = "BH")
```

<!-- badges: start -->

<!-- badges: end -->

Questions? Pls open an issue.
