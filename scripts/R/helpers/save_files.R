#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Save files ####

# Classification and metadata
write_csv(classification, "output/R/classification.csv")
metadata_clean <- metadata %>% 
  filter(Sample_ID %in% colnames(phage_tpm$contig))
write_csv(metadata_clean, "output/R/metadata_clean.csv")

## Metadata and classification variables
system("mkdir -p output/R/R_variables")
saveRDS(metadata, "output/R/R_variables/metadata.RDS")
saveRDS(classification, "output/R/R_variables/classification.RDS")
saveRDS(classification_gnmd, "output/R/R_variables/classification_gnmd.RDS")

## Taxon and core pies
system("mkdir -p output/R/taxon_pies")
for (taxlevel in names(core_pie)) {
  ggsave(paste0("output/R/taxon_pies/core_pie.", taxlevel, ".pdf"),
         core_pie[[taxlevel]], width = 12, height = 10)
}
for (tax in names(taxon_pie)) {
  ggsave(paste0("output/R/taxon_pies/taxon_pie.", tax,".pdf"), taxon_pie[[tax]], width = 12, height = 10)
}
for (tax in names(taxon_pie_overlap)) {
  ggsave(paste0("output/R/taxon_pies/taxon_pie_overlap.", tax,".pdf"), taxon_pie_overlap[[tax]], width = 8, height = 6)
}

write_lines(novel_families, "output/R/taxon_pies/novel_families.txt")
write_lines(novel_orders, "output/R/taxon_pies/novel_orders.txt")

for (core_or_all in names(classified_taxa)) {
  for (tl in names(classified_taxa[[core_or_all]])) {
    write_csv(classified_taxa[[core_or_all]][[tl]], paste0("output/R/taxon_pies/classified.", core_or_all, ".", tl, ".csv"))
  }
}

for (thing in names(pretty_pie$tibbles)) {
  for (tl in names(pretty_pie$tibbles[[thing]])) {
    wid <- 4.5
    hei <- 4.5
    if (tl == "Class") {
      wid <- 5.75
      hei <- 4.5
    }
    write_csv(pretty_pie$tibbles[[thing]][[tl]],
              paste0("output/R/taxon_pies/pretty_pie.", thing, ".", tl,".csv"))
    ggsave(paste0("output/R/taxon_pies/pretty_pie.", thing, ".", tl,".pdf"),
           pretty_pie$plots[[thing]][[tl]], width = wid, height = hei)
  }
}

for (tl in names(pretty_patches)) {
  ggsave(paste0("output/R/taxon_pies/pretty_patch.", tl, ".pdf"),
         pretty_patches[[tl]],
         width = 9, height = 4)
}

for (family_group in names(pretty_special_families)) {
  for (thing in names(pretty_special_families[[family_group]]$tibbles)) {
    write_csv(pretty_special_families[[family_group]]$tibbles[[thing]],
              paste0("output/R/taxon_pies/pretty_special_familes.", family_group, ".", thing, ".csv"))
    ggsave(paste0("output/R/taxon_pies/pretty_special_familes.", family_group, ".", thing, ".pdf"),
           pretty_special_families[[family_group]]$plots[[thing]],
           width = 3, height = 8)
  }
}

## aANI boxplot
ggsave("output/R/aANI_boxplot.pdf", aANI_boxplot, width = 6, height = 6)
write_csv(aANI_stats, "output/R/aANI_stats.csv")

## Completeness and genome length boxplots
system("mkdir -p output/R/completeness_and_genome_length")
for (metric in names(completeness_and_genome_length)) {
  ggsave(paste0("output/R/completeness_and_genome_length/", metric, ".pdf"), 
         completeness_and_genome_length[[metric]], width = 4, height = 4)
}
ggsave("output/R/completeness_and_genome_length/completeness_all_samples.pdf",
       completeness$plot, width = 5, height = 5)
write_csv(completeness$tibble, "output/R/completeness_and_genome_length/completeness_all_samples.csv")

## Alpha diversity plots and tables ####
system("mkdir -p output/R/alpha")
system("mkdir -p output/R/alpha/alpha_all")
system("mkdir -p output/R/alpha/alpha_core_or_not")
for (tax in names(alpha)) {
  ggsave(paste0("output/R/alpha/alpha_all/alpha.",tax,".pdf"),
         alpha[[tax]]$plot, width = 12, height=10)
  write_csv(alpha[[tax]]$table,
            paste0("output/R/alpha/alpha_all/alpha.",tax,".diversity.csv"))
  write_csv(alpha[[tax]]$kruskal,
            paste0("output/R/alpha/alpha_all/alpha.",tax,".kruskal.csv"))
  for (core_or_not in names(alpha_core_or_not)) {
    ggsave(paste0("output/R/alpha/alpha_core_or_not/alpha_core.",core_or_not,".",tax,".pdf"),
           alpha_core_or_not[[core_or_not]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_core_or_not[[core_or_not]][[tax]]$table,
              paste0("output/R/alpha/alpha_core_or_not/alpha_core.",core_or_not,".",tax,".diversity.csv"))
    write_csv(alpha_core_or_not[[core_or_not]][[tax]]$kruskal,
              paste0("output/R/alpha/alpha_core_or_not/alpha_core.",core_or_not,".",tax,".kruskal.csv"))
  }
}

# By country
system("mkdir -p output/R/alpha/alpha_by_country")
system("mkdir -p output/R/alpha/alpha_by_country_core_or_not")
for (country in names(alpha_by_country)) {
  for (tax in names(alpha_by_country[[country]])) {
    ggsave(paste0("output/R/alpha/alpha_by_country/alpha_by_country.",country,".",tax,".pdf"),
           alpha_by_country[[country]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_by_country[[country]][[tax]]$table,
              paste0("output/R/alpha/alpha_by_country/alpha_by_country.",country,".",tax,".csv"))
    for (core_or_not in names(alpha_core_or_not)) {
      ggsave(paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_by_country_core.",core_or_not,".",country,".",tax,".pdf"),
             alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$plot, width = 12, height=10)
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$table,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_by_country_core.",core_or_not,".",country,".",tax,".diversity.csv"))
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$kruskal,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_by_country_core.",core_or_not,".",country,".",tax,".kruskal.csv"))
    }
  }
}

# Absolute counts

ggsave("output/R/absolute_counts_all_samples.pdf", vlp_overview$plot,
       width = 8, height =4)
write_csv(vlp_overview$stats, "output/R/absolute_counts_all_samples.csv")

for (tax in names(alpha_abs)) {
  ggsave(paste0("output/R/alpha/alpha_all/alpha_abs.",tax,".pdf"),
         alpha_abs[[tax]]$plot, width = 12, height=10)
  write_csv(alpha_abs[[tax]]$table,
            paste0("output/R/alpha/alpha_all/alpha_abs.",tax,".csv"))
  for (core_or_not in names(alpha_abs_core_or_not)) {
    ggsave(paste0("output/R/alpha/alpha_core_or_not/alpha_abs_core.",core_or_not,".",tax,".pdf"),
           alpha_abs_core_or_not[[core_or_not]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_abs_core_or_not[[core_or_not]][[tax]]$table,
              paste0("output/R/alpha/alpha_core_or_not/alpha_abs_core.",core_or_not,".",tax,".diversity.csv"))
    write_csv(alpha_abs_core_or_not[[core_or_not]][[tax]]$kruskal,
              paste0("output/R/alpha/alpha_core_or_not/alpha_abs_core.",core_or_not,".",tax,".kruskal.csv"))
  }
}

# Absolute counts by country
for (country in names(alpha_abs_by_country)) {
  for (tax in names(alpha_abs_by_country[[country]])) {
    ggsave(paste0("output/R/alpha/alpha_by_country/alpha_abs_by_country.",country,".",tax,".pdf"),
           alpha_abs_by_country[[country]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_abs_by_country[[country]][[tax]]$table,
              paste0("output/R/alpha/alpha_by_country/alpha_abs_by_country.",country,".",tax,".csv"))
    for (core_or_not in names(alpha_abs_by_country_core_or_not)) {
      ggsave(paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_abs_by_country_core.",core_or_not,".",tax,".pdf"),
             alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$plot, width = 12, height=10)
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$table,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_abs_by_country_core.",core_or_not,".",tax,".diversity.csv"))
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$table,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_abs_by_country_core.",core_or_not,".",tax,".kruskal.csv"))
    }
  }
}

## Beta diversity plots and tables ####
system("mkdir -p output/R/beta/beta_all")
system("mkdir -p output/R/beta/beta_core_or_not")
for (tax in names(beta_plot_list)) {
  pcoa_path <- paste0("output/R/beta/beta_all/", tax, "_pcoa")
  pcoa_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_pcoa")
  hist_path <- paste0("output/R/beta/beta_all/", tax, "_hist")
  hist_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_hist")
  system(paste0("mkdir -p ", pcoa_path))
  system(paste0("mkdir -p ", hist_path))
  system(paste0("mkdir -p ", pcoa_path_core_or_not))
  system(paste0("mkdir -p ", hist_path_core_or_not))
  for (p in names(beta_plot_list[[tax]])) {
    if (p=="all") {
      wid <- 12
      hei <- 10
    } else {
      wid <- 14
      hei <- 4
    }
    for (q in names(beta_plot_list[[tax]][[p]])) {
      ggsave(paste0(pcoa_path,"/beta.", tax,".", p, ".", q, ".pcoa.pdf"),
             beta_plot_list[[tax]][[p]][[q]], width = wid, height = hei)
      for (core_or_not in names(beta_plot_list_core_or_not)) {
        ggsave(paste0(pcoa_path_core_or_not,"/beta_core.", core_or_not, ".", tax,".", p, ".", q, ".pcoa.pdf"),
               beta_plot_list_core_or_not[[core_or_not]][[tax]][[p]][[q]], width = wid, height = hei)
      }
      
      if (q == "control") { next }
      ggsave(paste0(hist_path,"/beta.", tax,".", p, ".", q, ".hist.pdf"),
             beta_dist[[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      for (core_or_not in names(beta_plot_list_core_or_not)) {
        ggsave(paste0(hist_path_core_or_not,"/beta_core.",core_or_not, ".", tax,".", p, ".", q, ".hist.pdf"),
               beta_dist_core_or_not[[core_or_not]][[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      }
    }
  }
  write_csv(rownames_to_column(beta_dist[[tax]]$dist_df, "Sample_ID"),
            paste0("output/R/beta/beta_all/beta_dist_", tax,".csv"))
  for (core_or_not in names(beta_plot_list_core_or_not)) {
    write_csv(rownames_to_column(beta_dist_core_or_not[[core_or_not]][[tax]]$dist_df, "Sample_ID"),
              paste0("output/R/beta/beta_core_or_not/beta_dist_core_or_not_", core_or_not, ".", tax,".csv"))
  }
}

# Absolute counts
for (tax in names(beta_abs_plot_list)) {
  pcoa_path <- paste0("output/R/beta/beta_all/", tax, "_pcoa")
  pcoa_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_pcoa")
  hist_path <- paste0("output/R/beta/beta_all/", tax, "_hist")
  hist_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_hist")
  system(paste0("mkdir -p ", pcoa_path))
  system(paste0("mkdir -p ", hist_path))
  system(paste0("mkdir -p ", pcoa_path_core_or_not))
  system(paste0("mkdir -p ", hist_path_core_or_not))
  for (p in names(beta_abs_plot_list[[tax]])) {
    if (p=="core") {
      width <- 17
      height <- 15
    } else {
      width <- 17
      height <- 5
    }
    for (q in names(beta_abs_plot_list[[tax]][[p]])) {
      ggsave(paste0(pcoa_path,"/beta_abs.", tax,".", p, ".", q, ".pcoa.pdf"),
             beta_abs_plot_list[[tax]][[p]][[q]], width = width, height = height)
      for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
        ggsave(paste0(pcoa_path_core_or_not,"/beta_abs_core.", core_or_not, ".", tax,".", p, ".", q, ".pcoa.pdf"),
               beta_abs_plot_list_core_or_not[[core_or_not]][[tax]][[p]][[q]], width = width, height = height)
      }
      
      if (q == "control") { next }
      ggsave(paste0(hist_path,"/beta_abs.", tax,".", p, ".", q, ".hist.pdf"),
             beta_abs_dist[[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
        ggsave(paste0(hist_path_core_or_not,"/beta_abs_core.",core_or_not, ".", tax,".", p, ".", q, ".hist.pdf"),
               beta_abs_dist_core_or_not[[core_or_not]][[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      }
    }
  }
  write_csv(rownames_to_column(beta_abs_dist[[tax]]$dist_df, "Sample_ID"),
            paste0("output/R/beta/beta_all/beta_abs_dist_", tax,".csv"))
  for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
    write_csv(rownames_to_column(beta_abs_dist_core_or_not[[core_or_not]][[tax]]$dist_df, "Sample_ID"),
              paste0("output/R/beta/beta_core_or_not/beta_abs_dist_core_or_not_", core_or_not, ".", tax,".csv"))
  }
}

# ## Heatmaps ####
system("mkdir -p output/R/heatmaps/contigs")
for (tax in names(order_heatmaps_row)) {
  ggsave(paste0("output/R/heatmaps/contigs/", tax, ".pdf"),
         order_heatmaps_row[[tax]],
         width = 20 + plotted_samples[[tax]]/5, height =2 + plotted_contigs[[tax]]/5,
         limitsize = FALSE)
}
system("mkdir -p output/R/heatmaps/tax_collapse")
for (tl in names(tax_collapse_heatmaps)) {
  system(paste0("mkdir -p output/R/heatmaps/tax_collapse/",tl))
  for (tax in names(tax_collapse_heatmaps[[tl]])) {
    tax_clean <- str_replace_all(tax, "\\|", "_")
    ggsave(paste0("output/R/heatmaps/tax_collapse/", tl,"/", tax_clean, ".pdf"),
           tax_collapse_heatmaps[[tl]][[tax]],
           width = 20 + plotted_samples_tax_collapse[[tl]][[tax]]/5, height =2+plotted_contigs_tax_collapse[[tl]][[tax]]/5,
           limitsize = FALSE)
  }
}

## TPM ####
system("mkdir -p output/R/relative_abundance/relative_abundance_overall/")
write_csv(phage_tpm$contig, "output/R/relative_abundance/phage_tpm.csv")

# Hostgroup
for (p in names(hostgroup_hist)) {
  ggsave(paste0("output/R/relative_abundance/relative_abundance_overall/hostgroup_hist_", p, ".pdf"),
         hostgroup_hist[[p]], width=8, height=4)
}
system("mkdir -p output/R/relative_abundance/relative_abundance_hostgroups/")
for (group in names(average_tpm_host_group)) {
  for (tl in names(average_tpm_host_group[[group]]$plots)) {
    wid <- 6
    hei <- 6
    if (tl == "Sample_ID") {
      wid <- 8
      hei <- 6
    }
    if (tl == "Country") {
      wid <- 20
      hei <- 40 
    }
    ggsave(paste0("output/R/relative_abundance/relative_abundance_hostgroups/average_TPM_Host_groups_core.", group, ".",tl,".pdf"),
           average_tpm_host_group[[group]]$plots[[tl]], width = wid, height = hei)
    write_csv(average_tpm_host_group[[group]]$tibbles[[tl]],
              paste0("output/R/relative_abundance/relative_abundance_hostgroups/average_TPM_Host_groups_core.", group, ".",tl,".csv"))
  }
}

# Core or not
for (core_or_not in names(core_or_not_hist)) {
  ggsave(paste0("output/R/relative_abundance/relative_abundance_overall/averageTPM_core_",core_or_not,".",".pdf"),
         core_or_not_hist[[core_or_not]], width=15, height=8)
}

## Average TPM
system("mkdir -p output/R/relative_abundance/relative_abundance_by_metavar/")
for (tax in names(average_tpm)) {
  for (metavar in names(average_tpm[[tax]]$plots)) {
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar/relative_abundance_", tax, ".", metavar, ".pdf"),
           average_tpm[[tax]]$plots[[metavar]],  width=40, height=15)
    write_csv(average_tpm[[tax]]$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_metavar/relative_abundance_", tax, ".", metavar, ".csv"))
  }
}

## Core TPM
system("mkdir -p output/R/relative_abundance/core_TPM/")
for (mvar in names(core_tpm_stats)) {
  wid <- 4
  hei <- 4
  if (mvar == "Country") {
    wid <- 6
    hei <- 6
  }
  ggsave(paste0("output/R/relative_abundance/core_TPM/core_TPM.", mvar, ".pdf"),
         core_tpm_plots[[mvar]], width = wid, height = hei)
  write_csv(core_tpm_stats[[mvar]], paste0("output/R/relative_abundance/core_TPM/core_TPM.", mvar, ".csv"))
}

system("mkdir -p output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/")
for (core_or_not in names(average_tpm_core_or_not_taxes)) { 
  for (tax in names(average_tpm_core_or_not_taxes[[core_or_not]])) {
    for (metavar in names(average_tpm_core_or_not_taxes[[core_or_not]][[tax]]$plots)) {
      hei <- 5
      wid <- 12
      if (metavar == "Sample_ID") {
        wid <- 90
      }
      ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/relative_abundance_", core_or_not, ".", tax, ".", metavar, ".pdf"),
             average_tpm_core_or_not_taxes[[core_or_not]][[tax]]$plots[[metavar]],  width = wid, height = hei, limitsize = FALSE)
      write_csv(average_tpm_core_or_not_taxes[[core_or_not]][[tax]]$tibbles[[metavar]],
                paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/relative_abundance_", core_or_not, ".", tax, ".", metavar, ".csv"))
    }
  }
}
for (thing in names(average_tpm_core_or_not)) {
  for (metavar in names(average_tpm_core_or_not[[thing]]$plots)) {
    wid <- 12
    hei <- 5
    if (metavar == "Sample_ID") {
      wid <- 20
      hei <- 40
    }
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/By_prevalence_", thing, "_relative_abundance.", metavar, ".pdf"),
           average_tpm_core_or_not[[thing]]$plots[[metavar]],  width = wid, height = hei, limitsize = FALSE)
    write_csv(average_tpm_core_or_not[[thing]]$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/By_prevalence_", thing,"_relative_abundance.", metavar, ".csv"))
  }
}

## Prevalence ####
system("mkdir -p output/R/prevalence/")
system("mkdir -p data/prevalence_tables/")
for (tl in names(prevalence_histo)) {
  if (tl=="Bee_pools") {
    wid=25
    plot <- prevalence_histo[[tl]]$plot # + scale_x_continuous(breaks = seq(3,150,3))
  } else if (tl=="Hives") {
    wid=12
    plot <- prevalence_histo[[tl]]$plot 
  } else {
    wid=8
    plot <- prevalence_histo[[tl]]$plot
  }
  ggsave(paste0("output/R/prevalence/prevalence.",tl,".pdf"),
         plot, width=wid, height=5)
  write_csv(prevalence_histo[[tl]]$table, paste0("output/R/prevalence/prevalence.",tl,".csv"))
  write_csv(presence_absence[[tl]], paste0("output/R/prevalence/presence_absence.",tl,".csv"))

  # write_csv(prevalence_histo[[tl]]$table, paste0("data/prevalence_tables/prevalence.",tl,".csv")) # This is also written to data to avoid backtracking. So it can be used at the top of the main script already.
}

## TPM per country ####
system("mkdir -p output/R/relative_abundance/relative_abundance_by_country")
for (core_or_not in names(average_tpm_per_country)) {
  for (tl in names(average_tpm_per_country[[core_or_not]])) {
    for (country in names(average_tpm_per_country[[core_or_not]][[tl]])) {
      ggsave(paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_core.", core_or_not, ".", tl, ".", country, ".pdf"),
             average_tpm_per_country[[core_or_not]][[tl]][[country]], height = 15, width = 35)
    }
  }
}
for (country in names(average_tpm_per_country_metavar)) {
  for (metavar in names(average_tpm_per_country_metavar[[country]]$plots)) {
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_metavar.", metavar, ".", country, ".pdf"),
           average_tpm_per_country_metavar[[country]]$plots[[metavar]], height = 5, width = 15)
    write_csv(average_tpm_per_country_metavar[[country]]$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_metavar.", metavar, ".", country, ".csv"))
  }
}

## Prevalence per country ####
system("mkdir -p output/R/prevalence/prevalance_by_country")
for (core_or_not in names(prevalence_plots_per_country)) {
  for (tl in names(prevalence_plots_per_country[[core_or_not]])) {
    for (country in names(prevalence_plots_per_country[[core_or_not]][[tl]])) {
      ggsave(paste0("output/R/prevalence/prevalance_by_country/prevalence_core.", core_or_not, ".", tl, ".", country, ".pdf"),
             prevalence_plots_per_country[[core_or_not]][[tl]][[country]],
             width = 30, height = 20)
    }
  }
}

## Upset ####
pdf(file="output/R/venns/country_upset.pdf",
    width = 25, height = 5) 
country_upset_plot
dev.off()


## Genera per family
system("mkdir -p output/R/gener_per_family/")
ggsave("output/R/genera_per_family/genera_per_family.pdf", genera_per_family_plot,
       width = 6, height = 6)
write_csv(genera_per_family_stats, "output/R/genera_per_family/genera_per_family.csv")
