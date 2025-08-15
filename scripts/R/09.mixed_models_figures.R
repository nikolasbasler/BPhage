library(tidyverse)
library(patchwork)
library(ggtext)

source("scripts/R/helpers/mixed_helpers.R")


#### Main figure
all_slopes <- read.delim("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm_vs_landuse.all_slopes.tsv") %>% 
  tibble() %>%
  rename(`Std. Error` = Std..Error)

plots <- list()
plots$gene_presence$Chitinase$`Fungicides and Bactericides` <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.Chitinase.Fungicides and Bactericides.rds")
plots$gene_presence$Chitinase$Herbicides <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.Chitinase.Herbicides.rds")
plots$gene_presence$`PAPS reductase`$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.PAPS reductase.Insecticides.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(0.76, 0.92),
                     labels = c(0.76, 0.92))
plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.DWV B.Cropland_in_2km_radius.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(24.84, 33.22),
                     labels = c(24.84, 33.22))
plots$pathogen_ct$BQCV$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.BQCV.Insecticides.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(23.37, 26.53),
                     labels = c(23.37, 26.53))
plots$gene_tpm$nosema_relabund <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/single_panel.RDS.Citinase.Vairimorpha_relabund.rds")
plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Cropland_in_2km_radius.rds")
plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Pesticides (total).rds")
plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Fungicides and Bactericides.rds")
plots$gene_tpm$`PAPS reductase`$Herbicides <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Herbicides.rds")

color_list <- list(dark = list(Chitinase = "#8B4513",
                               Glucosyltransferase = "#FFC300",
                               `PAPS reductase` = "#ef8f01",
                               PnuC = "#FFB89A",
                               Levanase = "#FF8C69",
                               `DWV B` = "#6A0DAD",
                               BQCV = "black"),
                   bright = list(Chitinase = "#8B4513",
                                 Glucosyltransferase = "#FFC300",
                                 `PAPS reductase` = "#ef8f01",
                                 PnuC = "#FFB89A",
                                 Levanase = "#FF8C69",
                                 `DWV B` = "#6A0DAD",
                                 BQCV = "black"))

common_legend_1 <- legend_factory(title = "Gene", 
                                  items = names(color_list$dark)[3],
                                  colors = unlist(color_list$dark)[3],
                                  position = "top")
common_legend_row2 <- legend_factory(title = "Gene",
                                     items = names(color_list$dark)[c(3,1)],
                                     colors = unlist(color_list$dark)[c(3,1)],
                                     position = "left")
common_legend_row3a <- legend_factory(title = "Pathogen",
                                      items = names(color_list$dark)[6:7],
                                      colors = unlist(color_list$dark)[6:7],
                                      position = "left")
common_legend_row3b <- legend_factory(title = "Gene",
                                      items = names(color_list$dark)[1],
                                      colors = unlist(color_list$dark)[1],
                                      position = "left") /
  plot_spacer() + plot_layout(heights = c(25,1))

paps_column <- plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides +
  plot_layout(ncol = 1, axes = "collect")

second_row <- 
  plots$gene_presence$`PAPS reductase`$Insecticides +
  plots$gene_presence$Chitinase$`Fungicides and Bactericides` +
  plots$gene_presence$Chitinase$Herbicides +
  common_legend_row2 +
  plot_spacer() +
  plot_layout(nrow = 1, axes = "collect", widths = c(2,2,2,1,1))

third_row <-   
  plots$gene_tpm$nosema_relabund +
  plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
  plots$pathogen_ct$BQCV$Insecticides +
  common_legend_row3b +
  common_legend_row3a +
  plot_spacer() +
  plot_layout(nrow = 1, axes = "collect", widths = c(3,3,3,1,1,1))


#### Supplementary figure

names_presence <- c(
  "RDS.Chitinase.Fung & Bact – Benzimidazoles.rds",
  "RDS.Chitinase.Fung & Bact – Dithiocarbamates.rds",
  "RDS.Chitinase.Fung & Bact - nes.rds",
  "RDS.Chitinase.Herbicides – Amides.rds",
  "RDS.Chitinase.Herbicides – Triazines.rds",
  "RDS.Chitinase.Herbicides – Urea derivates.rds",
  "RDS.Chitinase.Herbicides - nes.rds",
  "RDS.Chitinase.Insecticides – Carbamates.rds",
  "RDS.Chitinase.Insecticides – Pyrethroids.rds",
  "RDS.Glucosyltransferase.Insecticides – Organo-phosphates.rds",
  "RDS.PAPS reductase.Herbicides – Dinitroanilines.rds",
  "RDS.PAPS reductase.Insecticides - nes.rds"
)
names_tpm <- c(
  "RDS.Levanase.Herbicides – Amides.rds",
  "RDS.Levanase.Herbicides – Urea derivates.rds",
  "RDS.PAPS reductase.Fung & Bact – Dithiocarbamates.rds",
  "RDS.PAPS reductase.Fung & Bact - nes.rds",
  "RDS.PAPS reductase.Herbicides – Amides.rds",
  "RDS.PAPS reductase.Herbicides – Phenoxy hormone products.rds",
  "RDS.PAPS reductase.Herbicides – Triazines.rds",
  "RDS.PAPS reductase.Herbicides – Urea derivates.rds",
  "RDS.PAPS reductase.Herbicides - nes.rds",
  "RDS.PAPS reductase.Insecticides – Carbamates.rds",
  "RDS.PAPS reductase.Insecticides – Pyrethroids.rds",
  "RDS.PnuC.Insecticides.rds"
)
names_ct <- c(
  "RDS.BQCV.Insecticides - nes.rds",
  "RDS.DWV B.Herbicides – Dinitroanilines.rds"
)


supp_plots <- list()

for (file in names_presence) {
  supp_plots$presence[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/", file))
}
for (file in names_tpm) {
  supp_plots$tpm[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/", file))
}
for (file in names_ct) {
  supp_plots$ct[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/", file))
}

sup_figure_presence <- wrap_plots(supp_plots$presence, axes = "collect_y")
sup_presence_legend <- legend_factory(title = "Gene",
                                      items = names(color_list$dark)[1:3],
                                      colors = unlist(color_list$dark)[1:3],
                                      position = "top")

sup_figure_tpm <- wrap_plots(supp_plots$tpm, axes = "collect_y")
sup_tpm_legend <- legend_factory(title = "Gene",
                                 items = names(color_list$dark)[3:5],
                                 colors = unlist(color_list$dark)[3:5],
                                 position = "top")

sup_figure_ct <- wrap_plots(supp_plots$ct, axes = "collect_y")
sup_ct_legend <- legend_factory(title = "Pathogen",
                                items = names(color_list$dark)[6:7],
                                colors = unlist(color_list$dark)[6:7],
                                position = "top")

#### Save files

system("mkdir -p output/R/genes_pathogens_and_landuse/selected_graphs")

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/paps_tpm.pdf",
       paps_column, height = 12, width = 3.5)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/goi_presence.pdf",
       second_row, height = 3, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/nosema_relabund_and_ct.pdf",
       third_row, height = 3, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/paps_legend.pdf",
       common_legend_1, height = 1, width = 3)

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_presence.pdf",
       sup_figure_presence, height = 9, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_presence_legend.pdf",
       sup_presence_legend, height = 1, width = 5)

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_tpm.pdf",
       sup_figure_tpm, height = 9, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_tpm_legend.pdf",
       sup_tpm_legend, height = 1, width = 4)

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_ct.pdf",
       sup_figure_ct, height = 3, width = 6)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/sup_fig_ct_legend.pdf",
       sup_ct_legend, height = 1, width = 3)


