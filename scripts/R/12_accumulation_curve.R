library(vegan)
library(tidyverse)

metadata <- readRDS("data/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
# classification <- readRDS("output/R/R_variables/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()
present_in_all_countries <- read_lines("data/core_contigs.txt")

phages <- list()
phages$core <- present_in_all_countries
phages$all <- phage_tpm$contig


## Individual samples and different gut parts 
phage_presence_mat <- list()
gp_list <- list(mid = "mid", ile = "ile", rec = "rec", all = c("mid", "ile", "rec"))
for (set in names(phages)) {
  for (gut_set in names(gp_list)) {
    
    phage_presence_mat[[set]][[gut_set]] <- phage_tpm %>%
      pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
      left_join(., metadata[c("Sample_ID", "Gut_part")], by = "Sample_ID") %>%
      filter(
        Gut_part %in% gp_list[[gut_set]],
        contig %in% phages[[set]],
        # tpm != 0
      ) %>%
      group_by(contig, Sample_ID) %>%
      summarise(phage_presence = ifelse(sum(tpm) > 0, TRUE, FALSE), .groups = "drop") %>%
      pivot_wider(names_from = contig, values_from = phage_presence, values_fill = FALSE) %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix()
  }
}

## Bee pools
for (set in names(phages)) {
  phage_presence_mat[[set]]$Bee_pool <- phage_tpm %>%
    pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
    left_join(., metadata[c("Sample_ID", "Bee_pool")], by = "Sample_ID") %>%
    filter(
      contig %in% phages[[set]],
      # tpm != 0
    ) %>%
    group_by(contig, Bee_pool) %>%
    summarise(phage_presence = ifelse(sum(tpm) > 0, TRUE, FALSE), .groups = "drop") %>%
    pivot_wider(names_from = contig, values_from = phage_presence, values_fill = FALSE) %>%
    column_to_rownames("Bee_pool") %>%
    as.matrix()
  
}

## Accumulation curves
accumulation_plot <- list()
total_observed <- tibble()
for (set in names(phage_presence_mat)) {
  for (met_v in names(phage_presence_mat[[set]])) {
  
    accum <- specaccum(phage_presence_mat[[set]][[met_v]],
                       method       = "random",
                       permutations = 100)
    
    accum_tibble <- tibble(sites = accum$sites,
                           richness = accum$richness,
                           sd = accum$sd
    ) %>%
      mutate(lower = richness - sd,
             upper = richness + sd)

    total_observed <- tibble(slice = paste0(set, "_", met_v),
           total_observed = max(accum_tibble$richness)) %>%
      rbind(total_observed, .)
    
    accumulation_plot[[set]][[met_v]] <- ggplot(accum_tibble, aes(x = sites, y = richness)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_point(size = 0.5) +
      labs(
        x     = "Number of samples",
        y     = "Accumulated number of unique phages" ) +
      theme_minimal(base_size = 14) +
      ggtitle(paste0(set, " phages - gut part: ", met_v))
    
  }
}

system("mkdir -p output/R/accumulation_curves/")
for (set in names(accumulation_plot)) {
  for (met_v in names(accumulation_plot[[set]])) {
    ggsave(paste0("output/R/accumulation_curves/accumulation_curve.", set, ".", met_v,".pdf"),
           accumulation_plot[[set]][[met_v]],
           width = 6, height = 6)
  }
}
write_csv(total_observed, "output/R/accumulation_curves/total_observed.csv")



#
# 
# p <- accumulation_plot +
#   geom_smooth(
#     method  = "nls",
#     formula = y ~ SSasymp(x, Asym, R0, lrc),
#     se      = FALSE,      # turn off ggplot’s own CI (we already have lower/upper)
#     color   = "firebrick",
#     size    = 1
#   )
# 
# 
# smooth_layer_index <- length(p$layers)
# asym_fit <- p$layers[[smooth_layer_index]]$model
# class(asym_fit)

#
# 
# accum_df <- accum_tibble
# 
# N <- max(accum_df$sites)
# accum_df %>% filter(sites == N) %>% pull(richness)
# # Should be 2346 if you truly observe 2346 contigs in the pooled data.
# 
# # ─── 3. Fit a single Michaelis–Menten (asymp) curve via nls() ──────────────
# #    We'll use:     richness ~ Asym * sites / (Drop + sites)
# #    Good starting guesses:
# start_Asym <- max(accum_df$richness)      # e.g. ~2346
# start_Drop <- N * 0.2                     # e.g. ~20% of total samples
# 
# fit_mean_nls <- nls(
#   formula = richness ~ Asym * sites / (Drop + sites),
#   data    = accum_df,
#   start   = list(Asym = start_Asym, Drop = start_Drop),
#   control = nls.control(maxiter = 100)
# )
# 
# asymptote_est <- coef(fit_mean_nls)["Asym"]
# asymptote_est
# 
# grid_sites <- seq(1, N, length.out = 200)
# fitted_vals <- predict(fit_mean_nls, newdata = data.frame(sites = grid_sites))
# fit_df <- tibble(sites = grid_sites, fitted = fitted_vals)
# 
# ggplot(accum_df, aes(x = sites, y = richness)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.4) +
#   geom_point(size = 2, color = "steelblue") +
#   geom_line(size = 1, color = "steelblue") +
#   geom_line(data = fit_df, aes(x = sites, y = fitted),
#             color = "firebrick", size = 1) +
#   labs(
#     x        = "Number of samples",
#     y        = "Accumulated # of unique contigs",
#     title    = "Mean Accumulation Curve + Michaelis–Menten Fit",
#     subtitle = paste0("Estimated asymptote = ", round(asymptote_est, 1))
#   ) +
#   theme_minimal(base_size = 14)
# 
# 
# #
# 
# fit_asymp <- fitspecaccum(accum, model = "asymp")
# 
# asym_vec <- sapply(fit_asymp$models, function(m) {
#   # In case some permutations failed to converge, check for NULL or errors:
#   if (inherits(m, "nls")) {
#     coef(m)["Asym"]
#   } else {
#     NA
#   }
# })
# 
# mean_asym <- mean(asym_vec)
# sd_asym   <- sd(asym_vec)
# quantiles <- quantile(asym_vec, c(0.025, 0.5, 0.975))
# 
# ggplot(data.frame(Asym = asym_vec), aes(x = Asym)) +
#   geom_histogram(binwidth = 5, color = "black", fill = "salmon", alpha = 0.7) +
#   geom_vline(xintercept = mean_asym, color = "blue", linetype = "dashed") +
#   labs(
#     x = "Asymptote (total‐richness) from each permutation",
#     y = "Count",
#     title = "Distribution of ‘Asym’ across 100 specaccum permutations",
#     subtitle = paste0(
#       "Mean = ", round(mean_asym,1),
#       "; SD = ", round(sd_asym,1),
#       "; 95% CI [", round(quantiles[1],1), " – ", round(quantiles[3],1), "]"
#     )
#   ) +
#   theme_minimal()

