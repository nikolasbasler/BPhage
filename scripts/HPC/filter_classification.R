suppressMessages(library(readxl))
suppressMessages(library(tidyverse))

# Mind the order of the arguments!
# 1st argument: genomad's virus_summary.tsv
# 2nd argument: checkv's quality_summary.tsv
# 3rd argument: ICTV VMR excel sheet
# 4th argument: Output base name

args <- commandArgs(trailingOnly = TRUE)

# Load genomad output and append lowest classified taxon 
bphage_ALL_1kb_cross_95.85_virus_summary.tsv <- read.delim(args[1]) %>%
  mutate(lowest_taxon = str_split(taxonomy, ";") %>% sapply(function(x) tail(x, n = 1)))

# Load checkv output
quality_summary.tsv <- read.delim(args[2])

# Load ICTV taxonomy table
VMR_19_250422_MSL37 <- read_excel(args[3], 
                                  sheet = "VMR 19 MSL37", 
                                  col_types = c("numeric","numeric", "text", "text", "text",
                                                "text", "text", "text", "text", "text",
                                                "text", "text", "text", "text", "text",
                                                "text", "text", "text", "text", "text", 
                                                "text", "text", "text", "text", "text"))

# Make a vector of phage-exclusive taxa
#
# only_phages <- VMR_19_250422_MSL37 %>%
#   filter(`Host Source`=="bacteria")
# excluding_phages <- VMR_19_250422_MSL37 %>%
#   filter(!`Host Source`=="bacteria")
#
# Extending the search for bacterial and archeal hosts in my case only includes
# Caudoviricetes. No archaea virus order or family gets added.
only_phages <- VMR_19_250422_MSL37 %>%
  filter(`Host Source`=="bacteria" | `Host Source`=="archaea")
excluding_phages <- VMR_19_250422_MSL37 %>%
  filter(!`Host Source`=="bacteria" & !`Host Source`=="archaea")
phage_exclusive_taxa <- c()
for (level in c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family")) {
  taxa <- only_phages[level] %>%
    filter(!is.na(.data[[level]])) %>%
    unlist(use.names=FALSE) %>%
    unique() 
  for (tax in taxa) {
    if (!tax %in% unlist(excluding_phages[level])) {
      phage_exclusive_taxa <- c(phage_exclusive_taxa,tax)
    }
  }
}

# Filter genomad output for >= 50% complete phage genomes and write output csv
quality_summary.tsv %>%
  inner_join(., bphage_ALL_1kb_cross_95.85_virus_summary.tsv, by=join_by("contig_id"=="seq_name")) %>%
  filter(completeness>=50) %>%
  filter(lowest_taxon %in% phage_exclusive_taxa) %>%
  write_csv(paste0(args[4], "_phages.csv"))

# Filter genomad output for >= 50% complete genomes classified only as virus
# or "Unclassified" and write output csv
quality_summary.tsv %>%
  inner_join(., bphage_ALL_1kb_cross_95.85_virus_summary.tsv, by=join_by("contig_id"=="seq_name")) %>%
  filter(completeness>=50) %>%
  filter(lowest_taxon=="Unclassified" | lowest_taxon=="Viruses") %>%
  write_csv(paste0(args[4], "_unclassified_viruses.csv"))

# Filter genomad output for Picobirnaviridae irrespective of CheckV's completeness estimate
# (because they are segemented)
quality_summary.tsv %>%
  inner_join(., bphage_ALL_1kb_cross_95.85_virus_summary.tsv, by=join_by("contig_id"=="seq_name")) %>%
  filter(lowest_taxon=="Picobirnaviridae") %>%
  write_csv(paste0(args[4], "_picobirna.csv"))