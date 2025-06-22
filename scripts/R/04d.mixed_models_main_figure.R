library(tidyverse)
library(patchwork)

source("scripts/R/helpers/mixed_helpers.R")

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
# plots$gene_presence$`PAPS reductase`$`Herbicides – Dinitroanilines` <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.PAPS reductase.Herbicides – Dinitroanilines.rds")
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


#####
# PATCH 1

common_legend_genes <- legend_factory(title = "Gene", 
                                items = names(color_list$dark)[c(3,1)],
                                colors = unlist(color_list$dark[c(3,1)]),
                                position = "left")
common_legend_pathogens <- legend_factory(title = "Pathogen", 
                                      items = names(color_list$dark)[6:7],
                                      colors = unlist(color_list$dark[6:7]),
                                      position = "left")

common_legend_row1 <- legend_factory(title = "Gene", 
                                     items = names(color_list$dark)[3],
                                     colors = unlist(color_list$dark)[3],
                                     position = "left") /
  plot_spacer() + plot_layout(heights = c(50,1))
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


common_legend_row3 <- wrap_plots(plot_spacer(), common_legend_row3a, common_legend_row3b, plot_spacer(), nrow = 1)

first_row <- plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides +
  plot_layout(nrow = 1, axes = "collect")

paps_column <- plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides +
  plot_layout(ncol = 1, axes = "collect")

second_row <- 
  plots$gene_presence$`PAPS reductase`$Insecticides +
  plots$gene_presence$Chitinase$`Fungicides and Bactericides` +
  plots$gene_presence$Chitinase$Herbicides +
  # common_legend_row1 +
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


# main_figure <- common_legend_row1 / first_row /
#   common_legend_row2 / second_row /
#   third_row + plot_layout(heights = c(2,20,2,20,20))

main_figure <- first_row / second_row / third_row


#####
# PATCH 2


common_legend_1 <- legend_factory(title = "Gene", 
                                     items = names(color_list$dark)[3],
                                     colors = unlist(color_list$dark)[3],
                                     position = "top")
common_legend_2 <- wrap_plots(legend_factory(title = "Gene", 
                                     items = names(color_list$dark)[c(3,1)],
                                     colors = unlist(color_list$dark)[c(3,1)],
                                     position = "top"), plot_spacer(), widths = c(1,2))
common_legend_3a <- legend_factory(title = "Pathogen", 
                                      items = names(color_list$dark)[6:7],
                                      colors = unlist(color_list$dark)[6:7],
                                      position = "top")
common_legend_3b <- legend_factory(title = "Gene", 
                                      items = names(color_list$dark)[1],
                                      colors = unlist(color_list$dark)[1],
                                      position = "top")


# common_legend_row3 <- wrap_plots(common_legend_3b, common_legend_3a, nrow = 1)
common_legend_row3 <- wrap_plots(common_legend_3b, common_legend_3a, ncol = 1)

first <- 
  (plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides + plot_layout(axes = "collect")) /
  common_legend_1 +
  plot_layout(heights = c(100,1))

second <- 
  wrap_plots(plots$gene_presence$`PAPS reductase`$Insecticides +
  plots$gene_presence$Chitinase$`Fungicides and Bactericides` +
  plots$gene_presence$Chitinase$Herbicides +
  plot_layout(axes = "collect", nrow = 2),
  common_legend_2, ncol = 1, heights = c(100,1))

third <-   
  wrap_plots(
    plot_spacer() +
    plots$gene_tpm$nosema_relabund +
    plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
    plots$pathogen_ct$BQCV$Insecticides +
    plot_layout(nrow = 2, axes = "collect"),
    common_legend_row3, ncol = 1, heights = c(15,1))


system("mkdir -p output/R/genes_pathogens_and_landuse/selected_graphs")
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/wide.pdf",
       main_figure, height = 8, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/selected_graphs_3rd_row_legend.pdf",
       common_legend_row3, height = 3, width = 5)


ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/high1.pdf",
       first, height = 6, width = 6.5)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/high2.pdf",
       second, height = 6, width = 6.5)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/high3.pdf",
       third, height = 6.2, width = 6.5)

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/paps_col.pdf",
       paps_column, height = 12, width = 3.5)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/row2.pdf",
       second_row, height = 3, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/row3.pdf",
       third_row, height = 3, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/paps_legend.pdf",
       common_legend_1, height = 1, width = 3)


