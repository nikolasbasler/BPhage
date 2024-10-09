#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalenve Venn diagrams ####
# Output written to files already here.

venn_stats <- list()
tax_levels <- c("contig", "Genus", "Family")
met_v <- c("Season", "Gut_part", "Health")
for (tlvl in tax_levels) {
  for (m_var in met_v) {
    system(paste0("mkdir -p output/R/venns/core_yes/", tlvl, "/", m_var))
    system(paste0("mkdir -p output/R/venns/core_no/", tlvl, "/", m_var))
    system(paste0("mkdir -p output/R/venns/all/", tlvl, "/", m_var))
    for (cntr in levels(metadata$Country)) {
      venn_stats$all[[tlvl]][[m_var]][[cntr]] <- prevalence_venn(
        abtable = phage_ab_core_or_not$yes[[tlvl]],
        meta_var = m_var,
        country = cntr,
        subset = "all")
      for (core_or_not in names(phage_tpm_core_or_not)) {
        venn_stats[[core_or_not]][[tlvl]][[m_var]][[cntr]] <- prevalence_venn(
          abtable = phage_ab_core_or_not$yes[[tlvl]],
          meta_var = m_var,
          country = cntr,
          subset = paste0("core_",core_or_not))
      }
    }
  }
}