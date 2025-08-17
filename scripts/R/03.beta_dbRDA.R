library(RLdbRDA)
library(tidyverse)
library(patchwork)

source("scripts/R/custom_rldbrda.R")

metadata <- readRDS("data/metadata.RDS") %>%
  as.data.frame()
row.names(metadata) <- metadata$Sample_ID

core_or_not <- c("yes", "no")
plotting_lable <- list(yes = "core", no = "noncore", all = "all", 
                       abs_yes = "abs_core", abs_no = "abs_noncore", abs_all = "abs_all")

taxlevels <- c("contig", "Family")
beta_dists <- list()
beta_dists$Family$all <- read.csv("output/R/beta/beta_all/beta_dist_Family.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()
beta_dists$contig$all <- read.csv("output/R/beta/beta_all/beta_dist_contig.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()
beta_dists$Family$abs_all <- read.csv("output/R/beta/beta_all/beta_abs_dist_Family.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()
beta_dists$contig$abs_all <- read.csv("output/R/beta/beta_all/beta_abs_dist_contig.csv") %>%
  column_to_rownames("Sample_ID") %>%
  as.dist()

for (yes_no in core_or_not) {
  for (tax in taxlevels) {
    beta_dists[[tax]][[yes_no]] <- read.csv(paste0("output/R/beta/beta_core_or_not/beta_dist_core_or_not_", yes_no,".", tax,".csv")) %>%
      column_to_rownames("Sample_ID") %>%
      as.dist()
    beta_dists[[tax]][[paste0("abs_", yes_no)]] <- read.csv(paste0("output/R/beta/beta_core_or_not/beta_abs_dist_core_or_not_", yes_no,".", tax,".csv")) %>%
      column_to_rownames("Sample_ID") %>%
      as.dist()
  }
}

n <- 1
RDAs <- list()
RDA_plots <- list()
for (tax in taxlevels) {
  for (set in names(beta_dists[[tax]])) {
    print(paste0("Loop ", n, " of 12. Running RDA for tax: ", tax, ", and set: ", set, "."))
    meta_filt <- metadata %>%
      filter(Sample_ID %in% rownames(as.matrix(beta_dists[[tax]][[set]]))) %>% 
      select(Country, Season, Gut_part)
    
    if (str_detect(set, "abs_")) {
      meta_filt <- meta_filt %>%
        select(-Gut_part)
      custom_colors <- c("#8B4513", "#1C3A3A", "#FFA07A", "#1C3A3A")
    } else {
      custom_colors <- c("#ef8f01", "#1C3A3A", "#8B4513", "#1C3A3A", "#FFA07A", "#1C3A3A")
    }
    
    r_names_dist <- rownames(as.matrix(beta_dists[[tax]][[set]]))
    if (!identical(r_names_dist, rownames(meta_filt))) {
      print("ROWNAMES OF DIST MATRIX AND METADATA DONT MATCH! ABORTING")
      break
    }
    RDAs[[tax]][[set]] <- custom_rldbrda(beta_dists[[tax]][[set]], meta_filt)
    
    plot_data <- RLdbRDA::prepare_plot_data(RDAs[[tax]][[set]])
    
    p <- RLdbRDA::plot_dbrda(plot_data)
    p$data$bar_color <- custom_colors
    
    RDA_plots[[tax]][[set]] <- p + aes(fill = bar_color) + scale_fill_identity()
    n <- n + 1
  }
}

source("scripts/R/helpers/mixed_helpers.R")
common_legend <- legend_factory(title = "Thing", 
                                items = c("R2_A", "R2_B", "R2_C", "Cumulative R2"),
                                colors = c("#ef8f01", "#8B4513", "#FFA07A", "#1C3A3A"),
                                position = "bottom")

family_patch_horizontal <- 
  (RDA_plots$Family$all + RDA_plots$Family$no + RDA_plots$Family$yes) /
  common_legend +
  plot_layout(axes = "collect", heights = c(12,1))

family_patch_vertical <- RDA_plots$Family$all / RDA_plots$Family$no / RDA_plots$Family$yes / common_legend +
  plot_layout(axes = "collect", heights = c(rep(6,3),1))

# Save files
system("mkdir -p output/R/beta/beta_dbRDA/")

ggsave("output/R/beta/beta_dbRDA/dbRDA.Family_patch.horizontal.pdf", family_patch_horizontal, width = 10, height = 3)
ggsave("output/R/beta/beta_dbRDA/dbRDA.Family_patch.vertical.pdf", family_patch_vertical, width = 3.5, height = 9)
for (tax in taxlevels) {
  for (set in names(beta_dists[[tax]])) {
    write_delim(RDAs[[tax]][[set]], 
                paste0("output/R/beta/beta_dbRDA/dbRDA.", tax, ".", plotting_lable[[set]], ".tsv"),
                delim = "\t")
    ggsave(paste0("output/R/beta/beta_dbRDA/dbRDA.", tax, ".", plotting_lable[[set]], ".pdf"), RDA_plots[[tax]][[set]],
           width = 5, height = 5)
  }
}