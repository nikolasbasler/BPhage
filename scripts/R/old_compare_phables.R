library(tidyverse)


classification_gnmd <- read.csv("output/R/phage.filt.gnmd.classification.csv") %>%
  mutate(contig_length = contig_length/1000) %>%
  rename(length_kb = contig_length) %>%
  mutate(method = "tripple_ass")

number_of_contigs <- data.frame(method = "tripple_ass", contigs = nrow(classification_gnmd))

phage.filt.abundance.contig <- read.csv("output/R/phage.filt.abundance.contig.csv") %>%
  pivot_longer(-contig) %>%
  mutate(method = "tripple_ass")


phables_classification_gnmd <- read.csv("output/R/phables/phage.filt.gnmd.classification.csv") %>%
  mutate(contig_length = contig_length/1000) %>%
  rename(length_kb = contig_length)  %>%
  mutate(method = "phables")

phables_number_of_contigs <- data.frame(method = "phables", contigs = nrow(phables_classification_gnmd))

phables_phage.filt.abundance.contig <- read.csv("output/R/phables/phage.filt.abundance.contig.csv") %>%
  pivot_longer(-contig) %>%
  mutate(method = "phables")

contig_lengths <- rbind(classification_gnmd, phables_classification_gnmd) %>%
  ggplot(aes(x = method, y = length, fill = method)) +
  geom_boxplot() +
  ggtitle("Contig lengths") +
  theme(legend.position="none")
contig_lengths
completeness <- rbind(classification_gnmd, phables_classification_gnmd) %>%
  filter(Family != "Picobirnaviridae") %>%
  ggplot(aes(x = method, y = completeness, fill = method)) +
  geom_boxplot() +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggtitle("Completeness") +
  theme(legend.position="none")
completeness

num_contigs_and_mapped_reads <- rbind(phage.filt.abundance.contig, phables_phage.filt.abundance.contig) %>%
  group_by(method) %>%
  summarise(mapped_reads = sum(value)) %>%
  full_join(., rbind(number_of_contigs, phables_number_of_contigs), by = "method") %>%
  pivot_longer(-method) %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(name), scales = "free_y") +
  theme(legend.position="none")
num_contigs_and_mapped_reads

ggsave("output/R/phables/contig_lengths.pdf", contig_lengths,
       height = 6, width = 6)
ggsave("output/R/phables/completeness.pdf", completeness,
       height = 6, width = 6)
ggsave("output/R/phables/num_contigs_and_mapped_reads.pdf", num_contigs_and_mapped_reads,
       height = 6, width = 6)
  
