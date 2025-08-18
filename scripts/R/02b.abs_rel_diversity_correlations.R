library(tidyverse)
library(patchwork)

metadata <- readRDS("data/metadata.RDS") 

# Read data alpha diversity
alpha_family <- list()
alpha_family$all <- read.csv("output/R/alpha/alpha_all/alpha.Family.diversity.csv") %>% tibble()
alpha_family$noncore <- read.csv("output/R/alpha/alpha_core_or_not/alpha_core.no.Family.diversity.csv" )%>% tibble()
alpha_family$core <- read.csv("output/R/alpha/alpha_core_or_not/alpha_core.yes.Family.diversity.csv") %>% tibble()

# Count data alpha diversity
alpha_family_abs <- list()
alpha_family_abs$all <- read.csv("output/R/alpha/alpha_all/alpha_abs.Family.diversity.csv") %>% tibble()
alpha_family_abs$noncore <- read.csv("output/R/alpha/alpha_core_or_not/alpha_abs_core.no.Family.diversity.csv") %>% tibble()
alpha_family_abs$core <- read.csv("output/R/alpha/alpha_core_or_not/alpha_abs_core.yes.Family.diversity.csv") %>% tibble()

# Read data beta diversity
beta_family <- list()
beta_family$all <- read.csv("output/R/beta/beta_all/beta_dist_Family.csv") %>% tibble()
beta_family$noncore <- read.csv("output/R/beta/beta_core_or_not/beta_dist_core_or_not_no.Family.csv") %>% tibble()
beta_family$core <- read.csv("output/R/beta/beta_core_or_not/beta_dist_core_or_not_yes.Family.csv") %>% tibble()

# Count data beta diversity
beta_family_abs <- list()
beta_family_abs$all <- read.csv("output/R/beta/beta_all/beta_abs_dist_Family.csv") %>% tibble()
beta_family_abs$noncore <- read.csv("output/R/beta/beta_core_or_not/beta_abs_dist_core_or_not_no.Family.csv") %>% tibble()
beta_family_abs$core <- read.csv("output/R/beta/beta_core_or_not/beta_abs_dist_core_or_not_yes.Family.csv") %>% tibble()

# Shannon cor plotting function
correlate_shannons <- function(rel_alpha, abs_alpha, title) {
  shannon_relative <- rel_alpha %>%
    left_join(., metadata[c("Sample_ID", "VLPs_per_ul")], by = "Sample_ID") %>%
    filter(!is.na(VLPs_per_ul)) %>%
    rename(rel_Hill_Shannon = Hill_Shannon) %>%
    select(Sample_ID, rel_Hill_Shannon)
  
  shannon_absolute <- abs_alpha %>%
    rename(abs_Hill_Shannon = Hill_Shannon) %>%
    select(Sample_ID, abs_Hill_Shannon)
  
  p <- left_join(shannon_relative, shannon_absolute, by = "Sample_ID") %>%
    ggplot(aes(x = rel_Hill_Shannon, y = abs_Hill_Shannon)) +
    geom_point() +
    geom_smooth(formula = "y~x", method = "lm") +
    stat_cor(method = "spearman") +
    theme_minimal() +
    labs(x = "Hill Shannon (read data)",
         y = "Hill Shannon (count data)") #+
    # ggtitle(title)
  
  return(p)
}

# Bray cor plotting function
correlate_brays <- function(rel_beta, abs_beta, title) {
  bray_relative <- rel_beta %>%
    pivot_longer(-Sample_ID, names_to = "Sample_ID_2", values_to = "rel_Bray") %>%
    left_join(., metadata[c("Sample_ID", "VLPs_per_ul")], by = "Sample_ID") %>%
    filter(!is.na(VLPs_per_ul)) %>%
    select(-VLPs_per_ul) %>%
    left_join(., metadata[c("Sample_ID", "VLPs_per_ul")], by = join_by(Sample_ID_2 == Sample_ID)) %>%
    filter(!is.na(VLPs_per_ul)) %>%
    select(-VLPs_per_ul)
  
  bray_absolute <- abs_beta %>%
    pivot_longer(-Sample_ID, names_to = "Sample_ID_2", values_to = "abs_Bray")
  
  p <- left_join(bray_relative, bray_absolute, by = c("Sample_ID", "Sample_ID_2")) %>%
    ggplot(aes(x = rel_Bray, y = abs_Bray)) +
    geom_point(alpha = 0.1) +
    geom_smooth(formula = "y~x", method = "lm") +
    stat_cor(method = "spearman") +
    theme_minimal() +
    labs(x = "Bray-Curtis (read data)",
         y = "Bray-Curtus (count data)") # +
    # ggtitle(title)
  
  return(p)
}

titles <- c("all" = "all phages", "noncore" = "non-core", "core" = "core")

# Shannon correlation
shannon_cor_plots <- list()
for (set in c("all", "noncore", "core")) {
  shannon_cor_plots[[set]] <- correlate_shannons(
    rel_alpha = alpha_family[[set]],
    abs_alpha = alpha_family_abs[[set]],
    title = titles[[set]])
}

shannon_wrap <- wrap_plots(shannon_cor_plots, axes = "collect", ncol = 1)

# Bray correlation
bray_cor_plots <- list()
for (set in c("all", "noncore", "core")) {
  bray_cor_plots[[set]] <- correlate_brays(
    rel_beta = beta_family[[set]],
    abs_beta = beta_family_abs[[set]],
    title = titles[[set]]
  )
}

bray_wrap <- wrap_plots(bray_cor_plots, axes = "collect", ncol = 1)


# Save files

system("mkdir -p output/R/alpha/rel_abs_shannon_correlation")
system("mkdir -p output/R/beta/rel_abs_bray_correlation")

ggsave("output/R/alpha/rel_abs_shannon_correlation/alpha_rel_abs_shannon_correlation.wrap.pdf", 
       shannon_wrap, height = 8, width = 5)
ggsave("output/R/beta/rel_abs_bray_correlation/beta_rel_abs_bray_correlation.wrap.pdf", 
       bray_wrap, height = 8, width = 5)

for (set in c("all", "noncore", "core")) {
  ggsave(paste0("output/R/alpha/rel_abs_shannon_correlation/alpha_rel_abs_shannon_correlation.", set, ".pdf"),
         shannon_cor_plots[[set]]+ggtitle(set), height = 4, width = 4)
  ggsave(paste0("output/R/beta/rel_abs_bray_correlation/beta_rel_abs_bray_correlation.", set, ".pdf"),
         bray_cor_plots[[set]]+ggtitle(set), height = 4, width = 4)
}



