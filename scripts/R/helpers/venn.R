display_venn <- function(cat_list, f_name=NULL, ...){
  grid.newpage()
  venn_object <- venn.diagram(cat_list, filename = f_name, 
                              disable.logging=TRUE, 
                              ...)
}

find_intersects <- function(combos, taxes_list) {
  inter <- intersect(taxes_list[[combos[1]]], taxes_list[[combos[2]]])
  if (length(combos)>2) {
    for (l in 3:(length(combos))) {
      inter <- intersect(inter, taxes_list[[combos[l]]])
    }
  }
  return(inter)
}

prevalence_venn <- function(abtable, meta_var, country) {
  tax <- colnames(abtable)[1]
  
  if (country!="all") {
    abtable <- abtable %>%
      select(all_of(tax), starts_with(country))
  }
  ab_long_tibble <- abtable %>%
    pivot_longer(-all_of(tax)) %>%
    filter(value > 0) %>%
    left_join(., metadata, by=join_by(name==Sample_ID)) %>%
    select(all_of(c(tax, meta_var)))
  
  contained_taxes <- list()
  for (variable in levels(ab_long_tibble[[meta_var]])) {
    contained_taxes[[variable]] <- ab_long_tibble %>%
      filter(.data[[meta_var]]==variable) %>%
      filter(.data[[tax]]!="unclassified") %>%
      select(all_of(tax)) %>%
      unique() %>%
      unlist(use.names = FALSE)
  }
  n_cats <- length(levels(ab_long_tibble[[meta_var]]))
  colors <- c("#E69F00", "#56B4E9", "#009E73", "#999999")[1:n_cats]
  png_file <- paste0("output/R/venns/",tax,"/",meta_var,"/prevalence.Venn.",country,".",meta_var,".",tax,".png")
  
  display_venn(cat_list = contained_taxes, 
               fill = colors, 
               euler.d = FALSE, 
               scaled=FALSE,
               main = paste0(country," - ", meta_var," (",tax,")"),
               f_name = png_file
               )
  
  # Gather stats
  intersections <- list()
  in_higher_order_overlap <- NULL
  for (i in rev(2:n_cats)) {
    combs_matrix <- names(contained_taxes) %>%
      combn(i)
    for (j in 1:ncol(combs_matrix)) {
      comb <- combs_matrix[,j]
      key <- paste(comb, collapse=".")
      intersections[[key]] <- find_intersects(
        combos = comb,
        taxes_list = contained_taxes) %>%
        setdiff(., in_higher_order_overlap)
    }
    in_higher_order_overlap <- c(in_higher_order_overlap, intersections[[key]])
  }
  for (i in 1:length(contained_taxes)) {
    others <- contained_taxes[-i] %>%
      unlist(use.names = FALSE) %>%
      unique()
    key <- names(contained_taxes)[i]
    intersections[[key]] <- setdiff(contained_taxes[[i]], others)
  }
  
  stats_out <- NULL
  for (int in names(intersections)) {
    stats_out <- paste(intersections[[int]], collapse=", ") %>%
      paste(int, ., sep=": ") %>%
      c(stats_out, .)
  }
  stats_file <- paste0("output/R/venns/",tax,"/",meta_var,"/prevalence.Venn.",country,".",meta_var,".",tax,".txt")
  write_lines(stats_out, stats_file)
  
  return(intersections)
}
