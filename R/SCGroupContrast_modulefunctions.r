# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com


#' FitLmerContrast - for data with pre post treatment and 2 response groups; within each cell type contrast difference
#'  in fold change between groups, baseline difference, and fold change across groups of module scores.
#'
#' @param module_data_frame data for each cell to model -- designed to be scores for modules (columns) returned by scglmmr::WeightedModuleScore
#' @param lme4metadata metadata for model fit
#' @param f1 model furmula
#' @param contrast_list list of linear model contrasts.
#' @param plot_savepath path to save results
#' @param celltype_column the column in metadata with the cluster / celltype designation for each cell
#'
#' @return
#' @import ggplot2
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
#' @importFrom rlang sym
#' @importFrom emmeans emmeans
#' @importFrom ggpubr ggpaired rremove
#' @importFrom egg ggarrange
#' @importFrom broom tidy
#' @importFrom dplyr select mutate bind_rows if_else filter
#' @importFrom tidyr spread everything
#' @export
#'
#' @examples
#'\dontrun{
#' # load data from single cell data object
#' Seurat = readRDS("my_seurat_object.rds")
#'
#' # add cellwise module score for each signature
#' mod_scores = WeightedCellModuleScore(seurat_object = Seurat,
#'                                      module_list = btm,
#'                                      threshold = 0.1,
#'                                      return_weighted = FALSE, cellwise_scaling = FALSE,
#'                                      Seurat_version = "2")
#' Seurat = AddMetaData(Seurat,metadata = mod_scores)
#' module_n = names(sig_test)
#'
#' # set up module data frame
#' module_df = Seurat@meta.data %>% select(barcode_check, celltype_joint, module_n)
#'
#' # format metadata as factors group_id is order leveled for contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
#' md = Seurat@meta.data %>%
#'   mutate(group_id = factor(treat_time,  levels = c('pre_low', 'post_low', 'pre_high', 'post_high'))) %>%
#'   mutate(sampleid = factor(sampleid)) %>%
#'   select(barcode_check, celltype_joint, sampleid,  age, group_id)
#'
#' # Fit mixed model
#' plot_savepath = paste0(my_figure_save_path, "/marginalmeans/"); dir.create(plot_savepath)
#'
#' # specify any random intercept model e.g.
#' f1 = 'modulescore ~ age + group_id + (1|sampleid)'
#'
#' # fit sc mod mixed model on ewighted module scores.
#' mm_res = FitLmerContrast(module_data_frame = module_df,
#'                            celltype_column = 'celltype',
#'                              metadata = md,
#'                              fixed_effects = NULL,
#'                              lmer_formula = f1,
#'                              plotdatqc = TRUE,
#'                              figpath = 'your/file/path')
#'}
FitLmerContrast = function(module_data_frame, celltype_column = 'celltype', metadata,
                               fixed_effects = NULL, lmer_formula = NULL, plotdatqc = TRUE, figpath){

  # specify model formula
  if(!is.null(lmer_formula)){
    f1 = lmer_formula
  } else if(!is.null(fixed_effects)) {
    fixef = paste0(fixed_effects, collapse = " + ")
    f1 = paste("modulescore ~ 0 +", fixef, '+ group_id + (1|subjectid)')
  } else {
    f1 = paste("modulescore ~ 0 +", 'group_id + (1|subjectid)')
  }
  f_store = f1
  f1 = stats::as.formula(f1)

  # checks on minimal required metadata and formatting of group levels
  print('model specified:'); print(f1)
  stopifnot(is.factor(metadata$group_id))
  cat("predictor / perturbation variable must be a correctly ordered factor, testing: \n",
      "1: difference in treatmant effect between groups: ",
      levels(metadata$group_id)[3:4], " fold change  vs: ", levels(metadata$group_id)[1:2], "fold change  \n",
      "2: pre treatment baseline difference between groups: ",
      levels(metadata$group_id)[c(1,3)], " fold change \n",
      "3: treatment effect across both groups combined \n", "a prori contrasts across `ggroup_id` variable: "
  )
  # specify custom contrasts difference in treatment between groups,
  # treatment effect across groups, per-treatment difference between groups.
  c00 = c(1,0,0,0)
  c01 = c(0,1,0,0)
  c10 = c(0,0,1,0)
  c11 = c(0,0,0,1)
  contrast_list = list(
    "time1vs0_group2vs1" = (c11 - c10) - (c01 - c00),
    "time1vs0" = (c11 + c01) / 2 - (c10 + c00) / 2,
    "time0_group2vs1" = c10 - c00
  )
  print(contrast_list)

  # print data format note
  # print("module_data_frame columns: 1 = cell barcode 2 = celltype 3-nmodules = name of each module")
  #print("celltype column of lme4metadata and module_data_frame mmust be 'celltype_joint' ")

  # init storage
  res_celltype = res_list = list()

  # define model data
  module_names = names(module_data_frame)
  cts = as.character(unique(metadata[[celltype_column]]))

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))

    # subset meta
    metsub = metadata[metadata[[celltype_column]] == cts[u], ]
    meta_fit = colnames(metsub)
    stopifnot(c( 'group_id', 'subjectid') %in% meta_fit)

    # subset module dataframe to cells in celltype u from metadata above
    mod_df = module_data_frame[rownames(metsub), ]

    # merge module data and model metadata
    dat = cbind(metsub, mod_df)

    # indexing over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]; print(module)
      dat_fit = dat[ ,c(module, meta_fit)]
      colnames(dat_fit)[1] = "modulescore"

      # fit models
      gvar = rlang::sym('group_id')
      m1 = tryCatch(lme4::lmer(formula = f1, data = dat_fit), error = function(e)  return(NA))

      # calculate least squares means
      emm1 = tryCatch(
        emmeans::emmeans(object = m1, specs = ~ {{ group_id }}, data = dat_fit, lmer.df = "asymptotic"),
        error = function(e) return(NA)
        )

      # error checking on model convergence and least squares means
      if( suppressWarnings(is.na(m1))) {
        print(paste0(" could not fit  model for ", module_names[i]))
        res_list[[i]] = NA
      }
      if( suppressWarnings(is.na(emm1))) {
        print(paste0(module_names[i],  " check vcov matrix could not calculate marginal means  "))
        res_list[[i]] = NA
      } else{

        # format model results
        marginal_mean_df = broom::tidy(emm1)[ ,1:2] %>% tidyr::spread('group_id','estimate')
        names(marginal_mean_df) = paste0(names(marginal_mean_df), "_marginal_mean")

        # apply custom contrast
        contrast_fit = emmeans::contrast(emm1, method = contrast_list)
        contrastdf = broom::tidy(contrast_fit)
        contrast_names = names(contrast_list)
        contrast_colnames = sapply(contrast_names, function(x) paste0(colnames(contrastdf), x)) %>%
          as.data.frame() %>%
          unlist(use.names = F) %>%
          as.character()

        # single line for feature i ; to be rbound
        contrastdf = cbind(contrastdf[1, ], contrastdf[2, ], contrastdf[3, ])
        names(contrastdf) = contrast_colnames

        # flag singular fits
        singular_fit = m1@optinfo$conv$lme4$messages
        contrastdf = contrastdf %>%
          dplyr::mutate(singular_fit = dplyr::if_else(is.null(singular_fit), true = 0, false = 1))


        # add a Wilcoxon test of baseline difference between groups
        baseline_data = dat_fit %>% filter(group_id %in% levels(dat_fit$group_id)[c(1,3)])
        baseline_data$group_id = droplevels(baseline_data$group_id)
        wt = wilcox.test(modulescore ~ group_id, data = baseline_data) %>%
          broom::tidy() %>%
          rename(p.value.t0.wilcox = p.value)

        # merge all results and  format into a single line of a list to be r bound.
        contrastdf = cbind(contrastdf, wt)

        if(isTRUE(plotdatqc)){
          # reformat minimal marginal means
          p1 = plot(emm1) +
            coord_flip() +
            theme_classic() +
            theme(plot.title = element_blank()) +
            xlab("model estimated marginal means") +
            theme(axis.title.y = element_text(size = 7)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), axis.title.x = element_blank()) +
            scale_x_continuous(position = "left") +
            theme(axis.text.y = element_text(size = 4)) +
            theme(axis.line = element_line(colour = 'black', size = 0.2)) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


          ## add baseline and FC p value to plot
          t0p = formatC(contrastdf$p.value.t0.wilcox, format = "e", digits = 4)
          t1p = formatC(contrastdf$p.valuetime1vs0_group2vs1, format = "e", digits = 4)
          titleplot = paste0(cts[u]," marginal means ",module)

          p2 = ggplot(dat_fit, aes(x = group_id, y = modulescore)) +
            geom_violin(aes(fill = group_id), show.legend = FALSE) +
            stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", color = "black", size = 0.3) +
            scale_fill_manual(values = c("#1E90FFB3", "#191970B3", "#FF0000B3", "#CD0000B3")) +
            ggtitle(paste0(f_store, "\n", cts[u], " ", module, "\n",
                           "baseline wilcox p = ", t0p, "\n",
                           "lmer Fold Change Delta p = ", t1p)) +
            theme_bw(base_size = 11) +
            theme(axis.title.y =  element_text(size = 7)) +
            theme(plot.title = element_text(size = 7)) +
            ylab(paste0(module, "  score  " )) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7)) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), "cm"))

          # combine marginal means and data view
          p3 = egg::ggarrange(plots = list(p2,p1), nrow = 1, widths = c(3,1), draw = FALSE)
          ggsave(p3, filename = paste0(figpath,"VLN ", titleplot, ".pdf"), width = 3.5, height = 3.5)
        }
        # store result vector in result list i for celltype u
        res_list[[i]] =
          cbind(contrastdf, marginal_mean_df) %>%
          dplyr::mutate(module = module) %>%
          dplyr::mutate(formula = f_store) %>%
          dplyr::select(module, tidyr::everything())
      }
    }
    res_list = res_list[!is.na(res_list)]
    resdf = do.call(rbind, res_list)

    # store bound results for each module for celltype u in list element u
    res_celltype[[u]] = cbind(celltype = cts[u], resdf)
  }
  # bind results from each cell type into single dataframe
  resdf_full = do.call(rbind, res_celltype)
  return(resdf_full)
}


#' FitLmer - This is a simpler version of FitLmerContrast that makes no assumptions about structure of underlying
#'  data except that it can accomodate a mixed effects model formula and that there are multiple cell types to be
#'   separately fitted. fit a single cell mixed effects linear model. Designed to fit an aggregated gene module score.
#'   Fits a model to each cell cluster to test a grouping, treatment or combination factor; returns model fits
#'   for maximum flexibility
#'
#' @param module_data_frame data for each cell to model -- designed to be scores for modules (columns) returned by scglmmr::WeightedModuleScore
#' @param lme4metadata metadata for model fit
#' @param f1 model furmula
#' @param celltype_column the column in metadata with the cluster / celltype designation for each cell
#'
#' @return
#' @importFrom lme4 lmer
#' @importFrom stats as.formula
#' @importFrom tidyr spread everything
#' @export
#'
#' @examples
#'\dontrun{
# # s is a seurat object
#' s = SubsetData(s, ident.use = tc)
#' s = NormalizeData(s,normalization.method = 'LogNormalize',assay = 'RNA')
#'
#'
#' # format metadata as factors group_id is order leveled for:
#' # contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
#' md = s@meta.data %>%
#'   filter(timepoint %in% c(0,1)) %>%
#'   mutate(group_id = paste(group, timepoint, sep = '_')) %>%
#'   mutate(group_id = factor(group_id,  levels = c('0_0', '0_1', '1_0', '1_1'))) %>%
#'   mutate(subjectid = factor(sampleid)) %>%
#'   select(celltype, subjectid, age, group_id) %>%
#'   mutate(age = as.numeric(age)) %>%
#'   droplevels()
#'
#' # qc data to remove celltypes with no cells for some subhects at both timepoints
#' # keeps MLE more stable for the estimates of random intercept
#' ct.si = apply(table(md$celltype, md$subjectid) , 1, min)
#' c.keep = names(ct.si[ct.si > 7])
#' md = md[md$celltype %in% c.keep, ]
#'
#'
#' # add single cell weighted module scores
#' # split to standardize within cell type
#' ct.md = split(md, f = md$celltype)
#' mod_scores = lapply(ct.md, function(x){
#'   scglmmr::WeightedCellModuleScore(gene_matrix = s@assays$RNA@data[ ,rownames(x)],
#'                                    module_list = mtor.sigs,
#'                                    threshold = 0,
#'                                    # standardize within protein celltype
#'                                    cellwise_scaling = TRUE,
#'                                    return_weighted = FALSE )
#' })
#' ms = bind_rows(mod_scores)
#'
#' # correctly order rows after the split.
#' ms = ms[match(x = rownames(md), table = rownames(ms)), ]
#' stopifnot(all.equal(rownames(ms), rownames(md)))
#'
#'
#' # specify model
#' f1 = 'modulescore ~ group_id + age + (1|subjectid)'
#'
#' # fit sc mod mixed model on ewighted module scores.
#' mm_res.m1 = scglmmr::FitLmerContrast(module_data_frame = ms,
#'                                      celltype_column = 'celltype',
#'                                      metadata = md,
#'                                      lmer_formula = f1,
#'                                      plotdatqc = TRUE,
#'                                      fixed_effects = NULL,
#'                                      figpath = plot_savepath)
#'}
FitLmer = function(module_data_frame,
                   celltype_column = 'celltype',
                   metadata,
                   lme4.formula = NULL) {
  if (is.null(lme4.formula)) {
    stop('specify lme4.formula e.g. f <- ~ 0 + time.group + sex + age + (1|subjectid)')
  }
  if(Negate(celltype_column %in% colnames(metadata))){
    stop('celltype_column must be a variable column) in metadata')
  }

  f1 = stats::as.formula(lmer_formula)

  # init storage
  res_celltype = res_list = list()

  # define model data
  module_names = names(module_data_frame)
  cts = as.character(unique(metadata[[celltype_column]]))

  # indexed over cell type
  for (u in 1:length(cts)) {
    suppressMessages(gc())
    # celltype indexed
    print(paste0("fitting models for ", cts[u]))

    # subset meta
    metsub = metadata[metadata[[celltype_column]] == cts[u],]
    meta_fit = colnames(metsub)

    # subset module dataframe to cells in celltype u from metadata above
    mod_df = module_data_frame[rownames(metsub),]

    # merge module data and model metadata
    dat = cbind(metsub, mod_df)

    # indexing over modules
    for (i in 1:length(module_names)) {
      # module indexed ; select each module and fit model
      module = module_names[i]
      print(module)
      dat_fit = dat[, c(module, meta_fit)]
      colnames(dat_fit)[1] = "modulescore"

      # fit models
      m1 = tryCatch(
        lme4::lmer(formula = f1, data = dat_fit),
        error = function(e)
          return(NA)
      )

      # error checking on model convergence and least squares means
      if (suppressWarnings(is.na(m1))) {
        print(paste0(" could not fit  model for ", module_names[i]))
        res_list[[i]] = NA
      }
      # store result vector in result list i for celltype u
      res_list[[i]] = m1
    }

    res_list = res_list[!is.na(res_list)]
    # store bound results for each module for celltype u in list element u
    res_celltype[[u]] = list(res_list)
  }
  # bind results from each cell type into single dataframe
  names(res_celltype) = cts
  return(resdf_full)
}

