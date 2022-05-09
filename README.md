scglmmr *S*ample-level Single-*C*ell *G*eneralized *L*inear *M*ultilevel
*M*odels in *R*
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

An R package for implementing mixed effects models on single cell data
with complex experiment designs. The package is flexible and can
accomodate many experiment designs. It was developed for analysis of
multimodal single cell data from many individuals assayed pre and post
perturbation such as drug treatment, where each individual is nested
within one or more response groups. The methods herein allow one to
compare the difference in perturbation response effects between groups
while modeling variation in donor expression. It also has many wrappers
for downstream enrichment testing and visualization.

Please see vignettes.

## Installation

``` r
devtools::install_github(repo = "https://github.com/MattPM/scglmmr")
library(scglmmr)
```

<img src="man/figures/scglmmr.overview.png" />

**With this type of experiment design, we can’t just color umap plots
and try to find the effects.** We need statistical models.

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

An overview of methods provided:

### 1. pseudobulk aggregated models

These functions implement wrappers around `limma` for fitting fixed
effect linear models and the `dream` method from the `variancePartition`
package. The dream method is the only way to test differential
expression while accomodating ‘random’ or varying effects. This is
statistically necessary to account for non-independence when we have
perturbation experiments with repeated measurements from the same
donors. To enable linear models (e.g. modeling the mean with a normal
distribution) to be fit to gene counts, dream accounts for the mean
variance trend via incorporating `voom` observational weights.

### 2. single cell gene and gene module level mixed effects models

Test perturbation effect using a gene level Poisson mixed model.

### 3. single cell module level mixed effect models

Test perturbation effects and differences in perturbation responses
between groups at the gene module level.

### 4. Downstream enrichment testing and visualization

Wrappers around methods from the [fast set gene enrichment
(fgsea)](https://www.biorxiv.org/content/10.1101/060012v2#disqus_thread)
and
[clusterProfiler](https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
R packages.

### Philosophy

The scglmmr package considers each cluster/ cell type as a separate
‘experiment’ similar to fitting separate models to different FACS sorted
leukocyte subsets followed by RNAseq. Fitting models separately to each
subset provides maximum flexibility and avoids issues with e.g. modeling
mean variance trends or count distributions for cell type specific genes
in subsets that do not express the gene while still enabling comparison
of, for example, coherent perturbation effects for the same gene across
individuals between different cell clusters. This approach is
particularly well suited for CITE-seq data with cells clustered based on
normalized protein expression levels. Typically our workflow consists of
denoising ADT data using our method [dsb](https://github.com/niaid/dsb)
followed by modeling the group level perturbation response effects using
scglmmr.

**Experiment designs (within each cluster / celltype) supported by
scglmmr** Below is a 2 group repeated measures experiment. This data can
be accomodated by scglmmr. More simple experiment designs are also
supported, for example data with 2 groups but not repeated pre/post
treatment measurements.

| sample         |    sampleid    |      timepoint | Group          |      sex       |
|:---------------|:--------------:|---------------:|:---------------|:--------------:|
| 101_t0         |      101       |             d0 | good           |       F        |
| 101_t1         |      101       |             d1 | poor           |       F        |
| 102_t0         |      102       |             d0 | good           |       M        |
| 102_t1         |      102       |             d1 | poor           |       M        |
| … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) | … (x n donors) |

<!-- badges: start -->
<!-- badges: end -->

Questions? Pls open an issue.
