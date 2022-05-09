
# prepare built in package data

# hallmark
h = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
categ = unique(h$gs_name)
hallmark = sig_name = list()
for (i in 1:length(categ)) {
  hallmark[[i]] = h %>%
    as.data.frame() %>%
    dplyr::arrange(gs_name) %>%
    dplyr::filter(gs_name == categ[i]) %>%
    dplyr::select(gene_symbol) %>%
    purrr::set_names(categ[i])
  sig_name[i] = names(hallmark[[i]])
  hallmark[[i]] = hallmark[[i]][ ,1]
}
names(hallmark) = unlist(sig_name)

usethis::use_data(hallmark, compress = "xz")
