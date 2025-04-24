library(tidyverse)
library(patchwork)

source("scripts/R/helpers/mixed_helpers.R")

all_slopes <- read.delim("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm_vs_landuse.all_slopes.tsv") %>% 
  tibble() %>%
  rename(`Std. Error` = Std..Error)

plots <- list()
plots$gene_presence$Chitinase$`Fungicides and Bactericides` <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.Chitinase.Fungicides and Bactericides.rds")
plots$gene_presence$Chitinase$Herbicides <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.Chitinase.Herbicides.rds")
plots$gene_presence$`PAPS reductase`$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.PAPS reductase.Insecticides.rds")
plots$gene_presence$`PAPS reductase`$`Herbicides – Dinitroanilines` <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.PAPS reductase.Herbicides – Dinitroanilines.rds")
plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.DWV B.Cropland_in_2km_radius.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(24.84, 33.22),
                     labels = c(24.84, 33.22))
plots$pathogen_ct$BQCV$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.BQCV.Insecticides.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(23.37, 26.53),
                     labels = c(23.37, 26.53))
plots$gene_tpm$nosema_relabund <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/single_panel.RDS.Citinase.Varimorpha_relabund.rds")


#####
# Gene TPM plot rescaling

cropland_fraction <- read.csv("data/land_cover_results.csv") %>%
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)
# make FAOSTAT_added_data:
source("scripts/R/helpers/FAOstat_table.R")

cropland_and_FAO <- FAOSTAT_added_data %>%
  select(Country, Item, est_use_in_2k_radius) %>%
  pivot_wider(id_cols = Country, values_from = est_use_in_2k_radius, names_from = Item) %>%
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "Cropland_in_2km_radius")],. , by = "Country") %>%
  distinct() %>%
  select_if(~ !any(is.na(.)))

lowest_highest <- cropland_and_FAO %>%
  pivot_longer(-Country, names_to = "Item") %>%
  group_by(Item) %>%
  summarise(lowest = min(value),
            highest = max(value))

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


in_figure_gene_tpm <- c("PAPS reductase; Cropland_in_2km_radius", 
                        "PAPS reductase; Pesticides (total)",
                        "PAPS reductase; Fungicides and Bactericides", 
                        "PAPS reductase; Herbicides")

sig_tests <- all_slopes %>% 
  filter(test_name %in% in_figure_gene_tpm) %>%
  left_join(., lowest_highest, by = "Item") %>%
  mutate(effect = linear_effect_fun(s = Estimate, h = highest, l = lowest)) %>%
  group_by(gene) %>%
  mutate(y_stretching_factor = max(effect) / effect) %>%
  mutate(which_y_end_to_stretch = if_else(
    Estimate[y_stretching_factor == 1] > 0,
    "upper_end",
    "lower_end")) %>% 
  ungroup() %>%
  mutate(y_stretching_factor = case_when(
    which_y_end_to_stretch == "upper_end" & Estimate < 0 ~ 1/y_stretching_factor,
    which_y_end_to_stretch == "lower_end" & Estimate > 0 ~ 1/y_stretching_factor,
    .default = y_stretching_factor
  )) %>%
  select(-effect)

for (t_name in in_figure_gene_tpm) {
  
  log_tpm_test <- sig_tests %>%
    filter(test_name == t_name)
  
  tested_gene <- log_tpm_test$gene
  tested_item <- log_tpm_test$Item
  
  plots$gene_tpm[[tested_gene]][[tested_item]] <- 
    mixed_model_plot(filt_test_tibble = log_tpm_test,
                     transform_fun = linear_fun,
                     effect_fun = linear_effect_fun,
                     dark_col = color_list$dark[[tested_gene]],
                     bright_col = color_list$bright[[tested_gene]],
                     y_axis_label = "Log rel. gene abund.")
}


#####
# PATCH

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
                                     position = "top")
common_legend_row2 <- legend_factory(title = "Gene", 
                                     items = names(color_list$dark)[c(1,3)],
                                     colors = unlist(color_list$dark)[c(1,3)],
                                     position = "top")
common_legend_row3a <- legend_factory(title = "Pathogen", 
                                     items = names(color_list$dark)[6:7],
                                     colors = unlist(color_list$dark)[6:7],
                                     position = "left")
common_legend_row3b <- legend_factory(title = "Gene", 
                                      items = names(color_list$dark)[1],
                                      colors = unlist(color_list$dark)[1],
                                      position = "left")

# common_legend_row3 <- wrap_plots(common_legend_row3a, common_legend_row3b, plot_spacer(), nrow = 1, widths = c(2,1, 1))
common_legend_row3 <- wrap_plots(plot_spacer(), common_legend_row3a, common_legend_row3b, plot_spacer(), nrow = 1, widths = c(1, 4, 2, 1))
common_legend_row3 <- wrap_plots(plot_spacer(), common_legend_row3a, common_legend_row3b, plot_spacer(), nrow = 1)


first_row <- plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides +
  plot_layout(nrow = 1, axes = "collect")

second_row <- plots$gene_presence$Chitinase$`Fungicides and Bactericides` +
  plots$gene_presence$Chitinase$Herbicides +
  plots$gene_presence$`PAPS reductase`$Insecticides +
  plots$gene_presence$`PAPS reductase`$`Herbicides – Dinitroanilines` +
  plot_layout(nrow = 1, axes = "collect")

# third_row <- plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
#   plots$pathogen_ct$BQCV$Insecticides +
#   plots$gene_tpm$nosema_relabund +
#   plot_spacer() +
#   plot_layout(nrow = 1, axes = "collect")
#   # common_legend_genes +
#   # common_legend_pathogens +
#   # plot_layout(nrow = 1, axes = "collect", widths = c(7, 7, 7, 1, 3, 3))

third_row <- plot_spacer() +
  plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
  plots$pathogen_ct$BQCV$Insecticides +
  plots$gene_tpm$nosema_relabund +
  plot_spacer() +
  plot_layout(nrow = 1, axes = "collect", widths = c(1, 2, 2, 2, 1))

third_row <-   plot_spacer() +
  plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
  plots$pathogen_ct$BQCV$Insecticides +
  plots$gene_tpm$nosema_relabund +
  plot_spacer() +
  plot_layout(nrow = 1, axes = "collect", widths = c(1,20,20,20,19))

# common_legend_genes +
# common_legend_pathogens +
# plot_layout(nrow = 1, axes = "collect", widths = c(7, 7, 7, 1, 3, 3))

# main_figure <- first_row / plot_spacer() / second_row / plot_spacer()/ third_row + plot_layout(heights = c(20,1,20,1,20))

# main_figure <- common_legend_row1 / first_row / plot_spacer() / 
#   common_legend_row2 / second_row / plot_spacer() / 
#   common_legend_row3 / third_row + plot_layout(heights = c(2,20,1,2,20,1,2,20))

# main_figure <- common_legend_row1 / first_row /
#   common_legend_row2 / second_row /
#   common_legend_row3 / third_row + plot_layout(heights = c(2,20,2,20,2,20))

main_figure <- common_legend_row1 / first_row /
  common_legend_row2 / second_row /
  third_row + plot_layout(heights = c(2,20,2,20,20))

ggsave("output/R/genes_pathogens_and_landuse/selected_graphs.pdf",
       main_figure, height = 9, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs_3rd_row_legend.pdf",
       common_legend_row3, height = 3, width = 4)

