scglmmr *S*ample-level *S*ingle-cell *G*eneralized *L*inear *M*ultilevel
*M*odels in *R*
================
MPM

<!-- README.md is generated from README.Rmd. Please edit that file -->

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

# this needs to be rewritten

``` r
# # joint limma model on joint clusters H1 H5 
# suppressMessages(library(tidyverse))
# suppressMessages(library(Seurat))
# suppressMessages(library(magrittr))
# suppressMessages(library(variancePartition))
# suppressMessages(library(here))
# 
# # init cores 
# n_cores = 8
# 
# # get all DE functions
# source("de_workflow-master/full_differential_expression_workflow_functions.R")
# source("de_workflow-master/downstream_analysis_functions.r")
# 
# # dir = here("mid_res/2_H1_H5_joint_limmaDE/")
# datapath = here("your_file_path/generated_data/")
# figpath = here("your_file_path/figures/")
# dir.create(datapath); dir.create(figpath)
# 
# # read seurat object 
# Seurat = readRDS(file = "data/your_seurat_object.rds")
# 
# # DS model -- get celltypes ; remove low representation celltypes and doublets. 
# celltypes  = GetCelltypes(SeuratObject = Seurat, celltype = "celltype_joint")
# mytable = GetSubjectCelltypeTable(Seurat.Object = Seurat,celltype_column = "celltype_joint",sample_column = "sample")
# remove = mytable$`low representation celltypes`
# remove = c(remove, "DOUBLET")
# celltypes = celltypes[celltypes %ni% remove]
# 
# 
# # pool pseudobulk library raw counts 
# sl = MakePseudobulkList(seurat.object = Seurat, celltype_column = "celltype_joint", 
#                         sample_column = "sample",vector.of.celltypes = celltypes)
# saveRDS(sl, file = paste0(datapath,"raw_pseudobulk_list.rds"))
# 
# # make metadata and model matrix for contrasts fit. 
# d1 = MakeMetaTableFromSeurat(seurat_object = Seurat,sample_column = 'sample',variable_column = "cohort_timepoint",
#                              aggregate.data.list = sl, celltypes.vector = celltypes) 
# 
# 
# # create model metadata and contrast matrix 
# d1 = d1 %>% 
#   select(sample, cohort_timepoint) %>%
#   mutate(sampleid = str_sub(sample,-6,-4)) %>% 
#   droplevels() %>% 
#   column_to_rownames("sample")
# d1m = d1 %>% mutate_if(is.character, as.factor) %$% cohort_timepoint
# d1m = model.matrix(~0 + d1m)
# colnames(d1m) = str_sub(colnames(d1m), start = 4, end = -1)
# 
# # qc model rank and detect unused factor levels. 
# stopifnot(Matrix::rankMatrix(d1m) == ncol(d1m))
# stopifnot(any(colSums(d1m) == 0) == FALSE)
# 
# # make contrasts -- specify any arbitrary contrast 
# contrast_matrix  = makeContrasts(
#   # high vs low groups day 1 fold change 
#   t1_HighvsLow = (d1_high - d0_high) - (d1_low - d0_low),
#   ## timepoint only 
#   t1_vs_t0 = ((d1_high + d1_low)/2) - ((d0_high + d0_low)/2), 
#   ##low vs high
#   levels = colnames(d1m)
# )
# 
# # need this formatting for VariancePartition
# rownames(contrast_matrix) = paste0("cohort_timepoint", rownames(contrast_matrix))
# contrast_matrix = as.matrix(contrast_matrix) %>% as.data.frame()
# 
# # QC contrast matrix 
# variancePartition::plotContrasts(contrast_matrix)
# 
# # Filter matrix this edgeR function called within to get genes for each celltype https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# d1sum_dge = NormalizePseudobulk(aggregate.list =  sl, normalization.method = "RLE", model.metadataframe = d1,
#                                 comparison.group = "cohort_timepoint",minimum.gene.count = 5)
# 
# # run pseudobulk modelvoom weights with lme4 model with random intercept for subject, contrast day 1 fc in h5 vs h1. 
# cl = makeCluster(n_cores)
# doParallel::registerDoParallel(cl = cl)
# 
# # specify formula for lme4 / dream. 
# form = ~ 0 + cohort_timepoint + (1|sampleid)  
# 
# # optional; for dream lme4 pipeline, strip cell type from colnames of counts: 
# # d1sum_dge = lapply(d1sum_dge, function(x){ colnames(x$counts) = str_sub(colnames(x$counts), -6,-1) ; return(x)}) 
# 
# # Get voom observational level weights. 
# v1 = lapply(d1sum_dge, function(x){ x = voom(counts = x, design = d1m, normalize.method = "none", save.plot = TRUE, plot = FALSE)})
# 
# # fit model. 
# fit1 = lapply(v1, function(x){ dream(exprObj = x, formula = form, data = d1, L = contrast_matrix) })
# 
# 
# # save outputs. 
# names(d1sum_dge) = names(v1) = names(fit1) = celltypes 
# saveRDS(d1sum_dge, file = paste0(datapath,"/d1sum_dge_h1h5.rds"))
# model_list = list(v1, fit1)
# names(model_list) = c("voom", "fit")
# saveRDS(model_list, file = paste0(datapath,"/embedded_model_list.rds"))
# 
```

### Downstream analysis on pseudobulk results examples using provided convenience functions

``` r
### Extract contrast and ranke results 
## Note if using a version of variancepartition from bioconductor  can use versions of these functions in the more detailed section below 
# the "Raw" was added for a n earlier version of variancepartition that used ebayes on mixed model results .
# this is not valid so the "raw" versions of the functions for use with older Seurat versions calculate the ordinary t statistic manually 
##########
# Create results objects and run fGSEA on BTM. 
d1time_res = GetContrastResultsRaw(limma.fit.object.list = fit1,
                                   coefficient.number = 1, 
                                   contrast.name = "d1_h1vsh5_time", 
                                   celltypes.vector = celltypes)
d1time_rank = GetRankResultsRaw(contrast.result.raw.list = d1time_res) 
#########


##### Downstream analysis 
# run full gene set enrichment on gene ranks. 
gsea1 = RunFgseaOnRankList(rank.list.celltype = d1time_rank, pathways = btm)

# get tidy enrichment BTM data and make bubbleplot with positive and negative enrichmemt. 
score = RbindGseaResultList(gsea_result_list = gsea1,NES_filter = -Inf,padj_filter = 0.1)
GSEABubblePlot(rbind_gsea_result_dataframe = score,  include_negative = TRUE, 
               save_name = "/btm_0.1",
               save_path = figpath, 
               width = 7.9, height = 6.2
               ) 


# Run and plot hypergeometric test
# term_df = readRDS("signature_curation/cluster_profiler_termdf.rds") # this is included in the repo 
hyp_time = RunHypergeometricTest(result_list = d1time_res, TERM2GENE_dataframe = term_df,
                                 pval_threshold = 0.05, logFC_threshold = 0, usefdr_threshold = FALSE)
PlotHypergeometric(hyperg_result = hyp_time,p.adjust.filter = 0.05,genenumber_filter = 2, 
                   title = "hypergeometric test LI BTM, FDR 0.05, gene filter n > 2, genes with raw p value < 0.05", 
                   savepath = figpath,
                   savename = "d1time_hypergeometric_padj0.05.pdf",
                   height = 6, width = 8
                   )


## plot all DE genes from a given contrast in a heatmap. 
all_gene = do.call(rbind, d1time_res ) %$% gene %>% unique
mtx = GetGeneMatrix(result.list = d7time_res, 
                    stat_for_matrix = "logFC",
                    gene_subset = all_gene, 
                    pvalfilter = 0.01, 
                    logfcfilter = 0.25)
pheatmap::pheatmap(mtx,
                   color = viridis::viridis(n = 99, option = "B"),
                   breaks =seq(from = 0, to = 2,length.out = 99),
                   fontsize = 8, width = 4, height = 10, fontsize_row = 2.7, fontsize_col = 7, border_color = NA,  
                   filename = paste0(figpath,"d1_genes_heatmap.pdf"))


# plot only the leading edge genes from enrichment in a heatmap. 
le = GetLeadingEdgeFull(gsea.list = gsea1, padj.filter = 0.2,NES.filter = -Inf)
genesub = do.call(rbind, le) %$% gene %>% unique 
mtx2 = GetGeneMatrix(result.list = d7time_res, 
                    stat_for_matrix = "logFC",
                    gene_subset = genesub, 
                    pvalfilter = -Inf, 
                    logfcfilter = 0.2)
pheatmap::pheatmap(mtx2,
                   color = viridis::viridis(n = 99, option = "B"),
                   breaks =seq(from = 0, to = 2,length.out = 99),
                   fontsize = 8, width = 4, height = 10, fontsize_row = 2.7, fontsize_col = 7, border_color = NA,  
                   filename = paste0(figpath,"d1_LEgenes_heatmap.pdf"))
```

<!-- badges: start -->

<!-- badges: end -->

Questions? Please open an issue in the github.
