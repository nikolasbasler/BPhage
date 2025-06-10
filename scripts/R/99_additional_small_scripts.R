library(tidyverse)
library(patchwork)


##### 
## Calculating this factor of 25, regarding the value bee pollination services vs. the value of honey and beeswax:

# according to the federal reserve bank of st. louis (https://fred.stlouisfed.org/graph/?g=RJJP)
exchange_rate_2005 <- read.csv("data/Dollar_to_Euro_exchange_rates_2005_from_fred.stlouisfed.org.csv") %>%
  filter(!is.na(DEXUSEU)) %>%
  reframe(mean_rate = mean(DEXUSEU)) %>%
  unlist(use.names = FALSE)

pollination_value_euro <- 153e9/2 # Gallai et al 2009 estimate 153 billion euro by insect pollination in 2005. Reilly et al 2024 estimate honey bees contribute about half.
pollination_value_dollar <- pollination_value_euro * exchange_rate_2005

FAOSTAT_bee_product_value_en_5.30.2025 <- read.csv("data/FAOSTAT_bee_product_value_en_5-30-2025.csv") %>%
  tibble()

bee_product_value <- FAOSTAT_bee_product_value_en_5.30.2025 %>%
  select(Item, Value) %>%
  reframe(bee_product_value = sum(Value)*1000) %>%
  unlist(use.names = FALSE)

pollination_value_dollar / bee_product_value

###### 
## Bees in decline? WiP
FAOSTAT_data_en_5.30.2025 <- read.csv("~/Downloads/FAOSTAT_data_en_5-30-2025.csv") %>%
  tibble()

FAOSTAT_data_en_5.30.2025 %>%
  select(Area, Year, Value) %>%
  filter(!is.na(Value)) %>%
  group_by(Area) %>%
  mutate(total_per_area = sum(Value)) %>% 
  ungroup() %>%
  mutate(rank = dense_rank(total_per_area),
         reverted_rank = 1+abs(rank - max(rank))) %>%
  filter(Year >= 1990) %>%
  # filter(reverted_rank <= 10) %>%
  # ggplot(aes(x = Year, y = Value, color = Area)) +
  ggplot(aes(x = Year, y = Value)) +
  geom_line() +
  facet_wrap(~Area, scales = "free_y")



##### 
## MAPPING STATS SUMMARY:

pre_mapping_stats <- read.delim("output/bphage_viper_output/read_stats.tsv") %>%
  tibble() %>%
  rename("Sample_ID" = Sample) %>%
  filter(!str_detect(Sample_ID, "Blank")) %>%
  pivot_longer(-Sample_ID, names_to = "metric", values_to = "raw_reads") %>%
  mutate(metric = factor(metric, levels = c("Raw_R1_plus_R2", "Deduplicated_R1_plus_R2", 
                                            "Trimmed_R1_plus_R2", "Trimmed_unpaired", "Trimmed_total",
                                            "Hostout_R1_plus_R2", "Hostout_unpaired", "Hostout_total"))) %>%
  group_by(metric) %>%
  summarise(total = sum(raw_reads),
            average = mean(raw_reads),
            total_in_millions = total / 1000000,
            average_in_millions = average / 1000000) 

phage.filt.abundance.contig_long <- read.csv("output/R/phage.filt.abundance.contig.csv") %>%
  tibble() %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "filtered_mapped")
abundance.table_long <- read.csv("output/mapping_stats_phages/stats.phages.mapped_reads.csv") %>% 
  tibble() %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "unfilterd_mapped") %>%
  filter(!str_detect(Sample_ID, "Blank"))

post_mapping_stats <- full_join(phage.filt.abundance.contig_long, abundance.table_long, by = c("contig", "Sample_ID")) %>%
  group_by(Sample_ID) %>%
  summarise(unfilterd_mapped = sum(unfilterd_mapped),
            filtered_mapped = sum(filtered_mapped, na.rm = TRUE)) %>%
  reframe(total_unfiltered = sum(unfilterd_mapped),
          avaerage_unfiltered = mean(unfilterd_mapped),
          total_filtered = sum(filtered_mapped, na.rm = TRUE),
          average_filtered = mean(filtered_mapped, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "mapped_reads") %>%
  mutate(mapped_reads_in_millions = mapped_reads/1000000)

pre_mapping_stats
post_mapping_stats


#####
## BLASTn of phage genomes vs. INPHARED

alignment_plot <- function(blast_result) {
  blast_result %>%
    add_row(stitle="Query", sseqid="", qstart=1, qend = blast_result$qlen[1]) %>% # Add dummy row
    slice(c(n(), 1:(n() - 1))) %>% # Move dummy row to the bottom of the dataframe
    mutate(text = paste(stitle, sseqid)) %>%
    mutate(text = ifelse(nchar(text)>150, sub("^(.{0,})(.{145})$", "[...]\\2", text), text)) %>% # Shorten overly-long labels
    mutate(text = factor(text, levels = rev(unique(text)))) %>%
    ggplot(aes(x = qstart, xend = qend, y = text, yend = text, fill = pident)) +
    geom_segment(aes(color = pident), linewidth = 2) +
    labs(title = unique(blast_result$qseqid), x = "Position", y = NULL, color = "percent\nidentity\n") +
    scale_colour_gradient(limits = c(0, 100), na.value = "black")
}

present_in_all_countries <- read_lines("data/core_contigs.txt")
phages_blastn <- read.delim("output/phages_blastn.tsv") %>%
  tibble()

blast_of_core <- phages_blastn %>%
  filter(qseqid %in% present_in_all_countries)

plot_list <- list()
hits <- c()
longest_texts <- c()
tax_plot_wraps <- c()
for (cont in unique(blast_of_core$qseqid)) {
  filt_blast <- blast_of_core %>%
    filter(qseqid == cont)
  
  plot_list[[cont]] <-alignment_plot(blast_result = filt_blast)
  
  hits <- filt_blast$sseqid %>%
    unique() %>%
    length() %>%
    c(hits, .)
  longest_texts <- paste(filt_blast$stitle, filt_blast$stitle) %>%
    nchar() %>%
    max() %>%
    c(longest_texts, .)
}

wrap <- wrap_plots(plot_list) + plot_layout(ncol = 1, guides = "collect", heights=hits)

### Save wrap
ggsave("output/R/phage_blast.core.pdf", wrap,
       width = 7.5 + min(150, max(longest_texts))/15, height = length(plot_list) + sum(hits)/5,
       limitsize = FALSE
       )



