scglmmr *S*ample-level *S*ingle-cell *G*eneralized *L*inear *M*ultilevel
*M*odels in *R*
================
MPM

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

``` r
devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
library(scglmmr)
```

## Single cell within cluster perturbation response differential expression

The purpose of this workflow is to do within cluster differential
expression in a perturbation study where there are multiple measurements
made on the same donors. It can accomodate nested groups (e.g. poor vs
good outcome), 2 timepoints and can accomodate any fixed effects to
directly control for covariates.

Contact information:  
**general inquiry: John Tsang john.tsang AT nih.gov**  
**code for this project: Matt Mulè mulemp AT nih.gov or permanent email
mattmule AT gmail**  
For code related questions please open an issue in the issues tab

## Quickstart : 3 major modeling approaches

1.  single cell gene level poisson GLMM

2.  single cell module level GLMMs. Functions to calculate single cell
    weighted and non weighted module scores. Fit GLMM on these module
    scores with any lme4 random intercept model formula and by default
    test the group level differences at baseline and the differene in
    the fold change between the groups. Returns a plot of the least
    squares means of each group (residual mean of hte module score after
    accounting for covariates) with 95% CI.

3.  Convenience wrappers around the only pseudobulk DE method that
    accomodates mixed models and measurement uncertainty from pooling
    unqeual library sizes: limma + voom + lme4 implemented via the
    variancepartition package with the function dream().

Plus many downstream analysis wrappers to implement hypergeometric
testing with clusterprofiler, gene set enrichment with fGSEA, other
plotting functions etc.

## Method overview

### 1\. Single cell gene level DE testing

This does single cell level differential expression testing of genes by
fitting a Poisson mixed glm. You can have to specify a reduced and full
model formula to extract p values and effect sizes.

Recommended: only test a subset of genes within each celltype. Testing
30K genes will take a really long time and will artificially inflate FDR
corrected p values. Ideally you can select genes within each celltype by
running a variable gene detection method on each celltype / cluster,
saving the resulting genes in a list for each celltype.

The function accepts an argument ‘celltype\_genes\_test’ this is a named
list indexed by celltypes. the name of each element of the list must
match the celltype\_joint argument in metadata.

``` r
# Single cell mixed model. 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Seurat))
source("de_workflow-master/FitGLMM.R")

# set data path for save
datapath = here("your_save_path_relative_to_project_root")

# load subset of genes to test for each celltype 
genes_test = readRDS(file = "_subset_of_genes_per_celltype")
gene_union = genes_test %>% unlist %>% unique 

# read Seurat object (or any object)
s = readRDS("path_to_seurat_object") %>%
  SetAllIdent(id = "time_cohort") %>%
  SubsetData(ident.use = "d1", subset.raw = TRUE) 

# get metadata for day 1
md = s@meta.data %>%
  droplevels() %>%
  select(barcode_check, sampleid, timepoint,  celltype_joint, nUMI ) %>%
  mutate(sampleid = factor(sampleid)) %>%
  mutate(timepoint = factor(timepoint, levels = c("d0", "d1"))) %>% 
  rename(celltype = celltype_joint, barcode = barcode_check)

# # get gene data for day 1 cohort 
gene_data = Matrix::t(s@raw.data[gene_union, ])

## specify model formula; f2 = reduced model f1 = full model
f1 = 'gene ~ offset(log(nUMI)) + timepoint + (1|sampleid)'
f2 = 'gene ~ offset(log(nUMI)) + (1|sampleid)'
#f1 = 'modulescore ~ group_id + (1|sampleid)'

### test 
mmres = FitGLMM(gene_data = gene_data, lme4metadata = md, reduced_model = f2, full_model = f1, celltype_genes_test = genes_test)

### EXAMPLE output. 
#   celltype  gene    p.value     padj fixed_effect singular_fit fullmodel                                       reducedmodel                        
#   <chr>     <chr>     <dbl>    <dbl>        <dbl>        <dbl> <chr>                                           <chr>                               
# 1 BC_Naive  CD69   5.14e- 2 1.03e- 1       0.350             0 gene ~ offset(log(nUMI)) + timepoint + (1|samp… gene ~ offset(log(nUMI)) + (1|sampl…
# 2 BC_Naive  IRF1   3.89e- 4 1.56e- 3       0.863             0 gene ~ offset(log(nUMI)) + timepoint + (1|samp… gene ~ offset(log(nUMI)) + (1|sampl…
# 3 BC_Naive  STAT1  4.15e- 1 5.53e- 1       0.416             0 gene ~ offset(log(nUMI)) + timepoint + (1|samp… gene ~ offset(log(nUMI)) + (1|sampl…
# 4 CD14_Mono IRF1   2.10e-10 4.21e-10       0.982             0 gene ~ offset(log(nUMI)) + timepoint + (1|samp… gene ~ offset(log(nUMI)) + (1|sampl…
# 5 CD14_Mono STAT1  4.89e-10 4.89e-10       0.895             0 gene ~ offset(log(nUMI)) + timepoint + (1|samp… gene ~ offset(log(nUMI)) + (1|sampl…
```

### now you can run gene set enrichment on thee results, for example:

``` r
## enrichment on these results 

### plot enrichment with a hypergeometric test against the li BTM 
scd1 =  split( mmres , f = mmres$celltype )
scd1 = lapply(scd1, function(x){x %>%
    dplyr::select(gene, celltype, fixed_effect, padj)  %>%
    dplyr::rename('logFC'= fixed_effect,  'P.Value' = padj) })


# term_df = readRDS("signature_curation/cluster_profiler_termdf.rds")
hypsc = RunHypergeometricTest(result_list = scd1, TERM2GENE_dataframe = term_df,pval_threshold = 0.05, logFC_threshold = 0.1)  
PlotHypergeometric(hyperg_result = hypsc,
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
activity scores between groups within clusters. You can either test the
difference between the groups from a single timepoint, or if your data
has for example, pre and post treatment effects, you can test the
difference in fold changes in the groups, in other words this
accomodates *multiple timepoints, many donors with nested groups* using
mixed effects models and by default is set up to test the delta between
groups of the fold change in a module score at baseline vs post
perturbation using contrasts and the estimated least squares means after
correcting for covariates in the model with the emmeans package.

Here is a metadata example for a 2 group repeated measures experiment. A
random intercept model is required to account for autocorrelation in
each donors post perturbation response with their baseline level of
expression.

| sample  | sampleid | timepoint | Group | gender |
| :------ | :------: | --------: | :---- | :----: |
| 101\_t0 |   101    |        d0 | low   |   F    |
| 101\_t1 |   101    |        d1 | high  |   F    |
| 102\_t0 |   102    |        d0 | low   |   M    |
| 102\_t1 |   102    |        d1 | high  |   M    |

As a motivating toy example for implementing a mixed model on single
consider the simulated pathological case below exhibiting “Simpsons
Paradox” <https://en.wikipedia.org/wiki/Simpson%27s_paradox> (image from
wiki). Failing to account for the baseline expression unique to each
donor could result in model conclusing a gene’s post treatment effect is
the opposite direction\! (up vs down). In the example below, you can
think of the colored model as the “donor aware” model with random effect
for dubject ID. The black line would represent the model without the
random effect that incorrectly concludes the effect of the treatment
(slope) over time (x axis) on expression (y axis) is DOWN, when we see
that both donors INCREASED expression of the gene.
![image](https://user-images.githubusercontent.com/15280712/95596739-c4f94f80-0a1b-11eb-97bb-02bfafdf31bb.png)

#### Part I

You first fit a module score to each single cell. Options: *threshold*
what percent of genes in the module must a cell express to get a score
(i.e. of 0.2, if the cell has less than 20% of genes \>0 in the module
don’t score) *return\_weighted* whether to return the weighted module
score which is the average \* weight where weight = number of genes in
the cell in that module with non-zeroexpression *cellwise\_scaling*
should the scores be scaled across cells? Not recommended unless you are
testing only a single celltype. *Seurat\_version* if this is set to “2”
it will accomodate seurat v2, if set to “3” will accomodate v3

#### Part II fitting models:

*The FitModuleMixModel Function* To simultaneously test the difference
in the fold changes due to treatment (for example) the levels of the
combined group + time factor "group\_id have to be ordered as follows
(can be any names, these are just the order): group1 timepoint 0, group1
timepoint 1, group2 timepoint 0, group2 timepoint 1

It is also possible to use custom contrasts by changing the argument to
contrast\_list from contrast\_2 (automatically loaded) to something
else.

Specify any lme4 model formula in the argument to f1 f1 = f1 =
‘modulescore ~ age + gender + group\_id + (1|sampleid)’

This function will also automatically save a helpful plot of the least
squares means (the residual means of each group after accounting for
covariates in the model) of the module over the levels of each combined
group\_timepoint factor.

``` r

# load data
Seurat = readRDS("my_seurat_object.rds")

scglmmr::

# add cellwise module score for each signature 
mod_scores = WeightedCellModuleScore(seurat_object = Seurat,
                                     module_list = btm,
                                     threshold = 0.2,
                                     return_weighted = FALSE, cellwise_scaling = FALSE, 
                                     Seurat_version = "2") 
Seurat = AddMetaData(Seurat,metadata = mod_scores)
module_n = names(sig_test)

# set up module data frame 
module_df = Seurat@meta.data %>% select(barcode_check, celltype_joint, module_n) 

# format metadata as factors group_id is order leveled for contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
md = Seurat@meta.data %>% 
  mutate(group_id = factor(adjmfc.time,  levels = c('d0 low', 'd1 low', 'd0 high', 'd1 high'))) %>%
  mutate(sampleid = factor(sampleid)) %>% 
  select(barcode_check, celltype_joint, sampleid,  age, group_id)

# Fit mixed model 
plot_savepath = paste0(my_figure_save_path, "/marginalmeans/"); dir.create(plot_savepath)

# specify any random intercept model 
f1 = f1 = 'modulescore ~ age + group_id + (1|sampleid)'

# fit sc mod mixed model on ewighted module scores. 
mm_res = FitModuleMixModel(module_data_frame = module_df, lme4metadata = md, f1 = f1, contrast_list = contrast_2, plot_savepath = plot_savepath)

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
with the “models view” of the data in the estimated marginal means (in
the margin) and the “data” view of the actual module scores across each
level of
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

### 3\. “Pseudobulk” differential expression

aka limma + voom + lme4 + dream  
This provides convenience wrappers around variance partitions “dream”
method. It fits a weighted mixed linear mixed model implemented via the
dream function in the variancepartitionpackage which is a wrapper around
limma, voom, and lme4. It uses variance modeling at the observational
level (voom) to weight each observation (each gene sample combination)
which accounts for technical variation in the library sizes of
pseudobulk libraries derived from different numbers of cells per donor
within a cluster (cell type).

You can fit any lme4 random intercept formula, and fit a contrast to
extract effects of interest (see below)

This pipeline is inspired by 3 simulation studies. The first two
described statistical power (in particular, enhanced type I error
control) derived from pooling single cell mRNA reads into pseudobulk
libraries for differential expression testing. 10.1038/nmeth.4612
10.1101/713412

Based on these we adopted the limma+voom method for pseudobulk testing,
using the duplicatecorrelation function to account for donor-specific
baseline expression (analagous to a random effect, but using the genome
wide average covariance trend for each gene rather than a gene specific
random effect model). One observaiton on our data that motivated
adoption of limma-voom is the way in which voom observational weights
accounts for measurement uncertainty derived from unequal library sizes
(see Law et. al. genome biology - the original voom paper). In the case
of pseudobulk libraries derived from single cell RNAseq data,
interpolated model weights should correlate with library sizes within
each cell type. Indeed we confirmed that the log count per million and
square root standard deviation of genes had the expected monotonically
decreasing trend within each cell type and that the resulting
interpolated model weights across genes in a given sample were highly
correlated with both the number of cells used to create the library and
the sample’s total mRNA library size.

Later, a more recent report introduced the muscat R package and
simulation studies for testing within cluster differential expression.
10.1101/713412

This paper included simulation and validation studies of established
bulk RNAseq methods for within cluster differential expression analysis.
In that report, a pseudobulk method that performed particularly well was
using limma with voom observational weights, and mixed effects (i.e. the
ability to specify a random effect for donor baseline expression) with
lme4 model formula implemented through the “dream” method available
through the vairance partition package. We implemented a weighted mixed
effects model using this limma + voom + lme4 method within each cell
type using a random intercept for each donor with contrast coding to
test for group level fold change differences.

``` r
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(here))
'%ni%' = Negate('%in%')

#devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
library(scglmmr)

# load seurat or sce object etc. 
s = readRDS("path/seuratobject.rds")

# define counts and metadata and subset to cells above rm seurat object from workspace 
meta = s@meta.data
umi = s@assays$RNA@counts

# remove seurat object 
rm(s); gc()

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "lineage", sample_column = "sample")
tab$celltypes_remove

# remove cells prior to pseudobulk analysis 
cells_keep = meta %>% 
  rownames_to_column("bc") %>% 
  filter(celltype %ni% tab$celltypes_remove) %$% 
  bc

# subset data 
meta = meta[cells_keep, ]
umi = umi[ ,cells_keep]

# pseudobulk workflow 
pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample", celltype_col = "lineage", avg_or_sum = "sum")
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",variable_column = "cohort_timepoint", pseudobulklist = pb)
dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)

# custom a priori contrasts
# foldchange_difference =  difference in fold change between groups 
# time1_foldchange = overall time1 vs baseline fold change across both groups 
# baseline_groups - baseline difference between groups 
c_mat = makeContrasts(
  foldchange_difference = (cohort_timepoint1_1 - cohort_timepoint1_0) - (cohort_timepoint0_1 - cohort_timepoint0_0),
  time1_foldchange = (cohort_timepoint1_1 + cohort_timepoint0_1) / 2  - (cohort_timepoint1_0 + cohort_timepoint0_0) / 2, 
  baseline_groups = (cohort_timepoint1_0 - cohort_timepoint0_0),
  levels = colnames(designmat)
)
  
# fit mixed model for the multi timepoint contrasts 
fit = scglmmr::dreamMixedModel(dge_lists = dge, apriori_contrasts = TRUE, sample_column = 'sample',
                               cell_metadata = meta, contrast_matrix = c_mat, design_matrix = designmat, 
                               fixed_effects = c('age', 'gender', 'cohort_timepoint'), 
                               lme4_formula =  '~ 0 + age + gender + cohort_timepoint + (1|sampleid)', 
                               plotsavepath = figpath, version = "2", ncores = 4)

# fit simple linear model for the baseline group level contrast 
bl = scglmmr::RunVoomLimma(dgelists = dge, design_matrix = designmat, do_contrast_fit = T, my_contrast_matrix = c_mat[ ,3])
```

### Downstream analysis on pseudobulk results

example.

``` r
# hlmk = readRDS(file = here("signature_curation/hallmark.rds"))
figpath = "your/path"
test = scglmmr::GetRankResults(limma.fit.object.list = bl, coefficient.number = 1, "test")
res = scglmmr::GetContrastResults(limma.fit.object.list = bl, coefficient.number = 1, contrast.name = "test")

# note the name of these functions is different (to calculate the dinart t statistic) for models fit with dreamMixedModel
fit_res = scglmmr::GetContrastResultsRaw(limma.fit.object.list = fit, coefficient.number = 1,contrast.name = "foldchangedifference")
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


# make tidy average data 
av = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta,sample_col = "sample", celltype_col = "celltype", avg_or_sum = 'average')
le_expr = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av, gsea.list = hlmk_ctm0, padj.filter = 0.1,NES.filter = -Inf)


# example plot of sample level average leading edge genes annotated 
scglmmr::LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr,modulename = "MODULENAME",celltype_plot = "TCELL",
                      metadata = meta, metadata_annotate = c('group', 'timepoint', 'age', 'gender'),
                      sample_column = 'sample', returnmat = F, 
                      savepath = figpath, savename = "filename")

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



### Hypergeometric 
load(termdf)
hyp = scglmmr::RunHypergeometricTest(result_list = fit_res, TERM2GENE_dataframe = termdf,
                                     pval_threshold = 0.1,logFC_threshold = 0, usefdr_threshold = FALSE)
scglmmr::PlotHypergeometric(hyperg_result = hyp, p.adjust.filter = 0.1,genenumber_filter = 2,
                            savepath = figpath,savename = "name", title = "title")

# calculate average module z score across samples (z score across samples of gene z scores)
# baseline_samples = c('')
lapply(av, function(x), x[ ,baseline_samples])
mz = scglmmr::AverageSampleModuleZscore(average.metacell.list = av,module.list = hlmk,use.module.subset = F)
```

<!-- badges: start -->

<!-- badges: end -->

Questions? Pls open an issue.
