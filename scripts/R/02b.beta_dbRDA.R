library(RLdbRDA)
source("scripts/R/custom_rldbrda.R")
library(tidyverse)
library(patchwork)

metadata <- readRDS("output/R/R_variables/metadata.RDS")

core_or_not <- c("yes", "no")
plotting_lable <- list(yes = "core", no = "noncore", all = "all")

taxlevels <- c("contig", "Family")
beta_dists <- list()
beta_dists$Family$all <- read.csv("output/R/beta/beta_all/beta_dist_Family.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()
beta_dists$contig$all <- read.csv("output/R/beta/beta_all/beta_dist_contig.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()
for (yes_no in core_or_not) {
  for (tax in taxlevels) {
    beta_dists[[tax]][[yes_no]] <- read.csv(paste0("output/R/beta/beta_core_or_not/beta_dist_core_or_not_", yes_no,".", tax,".csv")) %>%
      column_to_rownames("Sample_ID") %>%
      as.dist()
  }
}

RDAs <- list()
RDA_plots <- list()
for (tax in taxlevels) {
  for (set in names(beta_dists[[tax]])) {
    meta_filt <- metadata %>%
      filter(Sample_ID %in% rownames(as.matrix(beta_dists[[tax]][[set]]))) %>% 
      select(Country, Season, Gut_part)
    r_names_dist <- rownames(as.matrix(beta_dists[[tax]][[set]]))
    if (!identical(r_names_dist, rownames(meta_filt))) {
      print("ROWNAMES OF DIST MATRIX AND METADATA DONT MATCH! ABORTING")
      break
    }
    RDAs[[tax]][[set]] <- custom_rldbrda(beta_dists[[tax]][[set]], meta_filt)
    
    plot_data <- RLdbRDA::prepare_plot_data(RDAs[[tax]][[set]])
    
    RDA_plots[[tax]][[set]] <- RLdbRDA::plot_dbrda(plot_data) + 
      scale_fill_manual(values=c("#ef8f01", "#8B4513"),
                        labels=c(bquote(R^2), bquote('Cumulative' ~ R^2))) +
      ggtitle(paste0(tax, "_", plotting_lable[[set]])) +
      theme(legend.position = "bottom")
  }
}

family_patch_horizontal <- RDA_plots$Family$all + RDA_plots$Family$no + RDA_plots$Family$yes +
  plot_layout(guides = 'collect', axes = "collect") & theme(legend.position = "bottom")

family_patch_vertical <- RDA_plots$Family$all / RDA_plots$Family$no / RDA_plots$Family$yes +
  plot_layout(guides = 'collect', axes = "collect") & theme(legend.position = "bottom")

# Save files
system("mkdir -p output/R/beta/beta_dbRDA/")

ggsave("output/R/beta/beta_dbRDA/dbRDA.Family_patch.horizontal.pdf", family_patch_horizontal, width = 13, height = 5)
ggsave("output/R/beta/beta_dbRDA/dbRDA.Family_patch.vertical.pdf", family_patch_vertical, width = 3.5, height = 9.5)
for (tax in taxlevels) {
  for (set in names(beta_dists[[tax]])) {
    write_delim(RDAs[[tax]][[set]], paste0("output/R/beta/beta_dbRDA/dbRDA.", tax, ".", plotting_lable[[set]], ".tsv"))
    ggsave(paste0("output/R/beta/beta_dbRDA/dbRDA.", tax, ".", plotting_lable[[set]], ".pdf"), RDA_plots[[tax]][[set]],
           width = 5, height = 5)
  }
}

