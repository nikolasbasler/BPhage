library(tidyverse)


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



