## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")
## created with project : thy
# pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sample", celltype_col = "lineage", avg_or_sum = "sum")
# designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sample",variable_column = "cohort_timepoint", pseudobulklist = pb)
# dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat, minimum.gene.count = 5)
#
# # custom a priori contrasts time effect for the IRAE group vs non irae group and overall time effect and baseline between groups.
# c_mat = makeContrasts(
#   irae_delta = (cohort_timepoint1_1 - cohort_timepoint1_0) - (cohort_timepoint0_1 - cohort_timepoint0_0),
#   time1_delta = (cohort_timepoint1_1 + cohort_timepoint0_1) / 2  - (cohort_timepoint1_0 + cohort_timepoint0_0) / 2,
#   baseline_irae = (cohort_timepoint1_0 - cohort_timepoint0_0),
#   levels = colnames(designmat)
# )
#
# # fit mixed model
# fit = scglmmr::dreamMixedModel(dge_lists = dge, apriori_contrasts = TRUE, sample_column = 'sample',
#                                cell_metadata = meta, contrast_matrix = c_mat, design_matrix = designmat,
#                                fixed_effects = c('cohort_timepoint'),
#                                lme4_formula =  '~ 0 + cohort_timepoint + (1|sampleid)',
#                                plotsavepath = figpath, version = "2", ncores = 4)
# testdat = lapply(fit, function(x) x[1:40, ])
# testpb = list()
# for (i in 1:length(pb)) {
#   testpb = pb[[i]][rownames(testdat)[1:40], ]
# }
# saved testpb and testdat to scglmmr root

# testfit = readRDS("test_fit.rds")
# testpb = readRDS("test_pb.rds")
usethis::use_data(testfit, compress = "xz")
usethis::use_data(testpb, compress = "xz")




