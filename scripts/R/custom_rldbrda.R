# NMDS scree plot
NMDS_scree <- function(x) { # where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, vegan::metaMDS(x, autotransform = F, k = 1, trace = F)$stress),
    xlim = c(1, 10), ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress",
    main = "NMDS stress plot"
  )
  for (i in 1:10) {
    points(rep(i + 1, 10), replicate(10, vegan::metaMDS(x, autotransform = F, k = i + 1, trace = F)$stress))
  }
}

#' Functions are adapted or copied (because they were not exported) from https://github.com/raeslab/RLdbRDA/tree/main
custom_rldbrda <- function(distmat, meta, p_cutoff = 0.05) {
  r2 <- get_r2(distmat, meta)

  sign_r2 <- rownames(r2[which(r2$padj < p_cutoff), ]) # selects only variables significant in step 1

  if (length(sign_r2) < 1) {
    message("No significant features found!")
    return(NULL)
  }

  cumul <- get_cumul(distmat, meta[, sign_r2, drop = F])

  out <- combine_data(r2, cumul)

  return(out)
}

get_r2_single <- function(distmat, meta, feature) {
  capsc <- vegan::capscale(distmat ~ meta[, feature], na.action = na.omit)
  an <- vegan::anova.cca(capsc, permutations = 9999)

  Fa <- an["F"][[1]][[1]]
  r2 <- vegan::RsquareAdj(capsc)[[1]]
  r2adj <- vegan::RsquareAdj(capsc)[[2]]
  N <- nrow(na.exclude(meta[, feature, drop = FALSE]))
  pval <- an["Pr(>F)"][[1]][[1]]

  output <- cbind(feature, Fa, r2, r2adj, N, pval)

  return(output)
}

get_r2 <- function(distmat, meta) {
  features <- colnames(meta)

  all <- c()

  for (feature in features) {
    es <- get_r2_single(distmat, meta, feature)
    all <- rbind(all, es)
  }

  all <- data.frame(all)
  all$padj <- p.adjust(all$pval, method = "BH")

  rownames(all) <- all$feature

  return(all)
}

get_cumul <- function(distmat, meta) {
  mod0 <- vegan::capscale(distmat ~ 1) # H0: unconstrained ordination
  mod1 <- vegan::capscale(distmat ~ ., data = meta) # H1: full constrained ordination, all metadata


  attach(meta)

  step.res <- vegan::ordiR2step(mod0, scope = formula(mod1), data = meta, direction = "forward", Pin = 1, R2scope = TRUE, pstep = 100, permutations = 9999, trace = F) # forward stepwise dbRDA
  res <- step.res$anova


  row.names(res) <- gsub(pattern = "\\+ ", "", row.names(res))
  colnames(res) <- gsub(pattern = "Pr\\(>F\\)", "pval", colnames(res)) # replace column name
  colnames(res) <- paste0("RDAcumul_", colnames(res))
  res[, "RDAcumul_N"] <- nrow(meta)

  detach(meta)


  return(res)
}

combine_data <- function(r2, cumul) {
  all <- data.frame(merge(r2, cumul, by = "row.names", all = T), row.names = 1)
  all <- all[order(all$r2, decreasing = TRUE), ]
  all <- all[order(all$RDAcumul_R2.adj), ]

  return(all)
}

prepare_plot_data <- function(dbrda_data) {
  # Define the order of the y-axis
  plot_order <- rev(row.names(dbrda_data))

  # Impute missing values in RDAcumul_R2.adj with the largest non-NA value
  insignificant_features <- row.names(dbrda_data)[is.na(dbrda_data$RDAcumul_R2.adj)]
  max_val <- max(dbrda_data$RDAcumul_R2.adj, na.rm = TRUE)
  dbrda_data$RDAcumul_R2.adj[is.na(dbrda_data$RDAcumul_R2.adj)] <- max_val

  # Filter the data and reshape it to long form
  plot_data <- dbrda_data %>%
    rownames_to_column(var = "rowname") %>%
    mutate(
      r2adj = as.numeric(r2adj),
      rowname2 = str_replace(rowname, "[0-9]", ""),
      RDAcumul_R2.adj = ifelse(padj <= 0.05, RDAcumul_R2.adj, NA)
    ) %>%
    filter(rowname2 != "<All variables>") %>%
    select(r2adj, RDAcumul_R2.adj, Stage, rowname2, rowname) %>%
    pivot_longer(c(-rowname, -rowname2, -Stage), names_to = "variable", values_to = "value") %>%
    mutate(significant = ifelse(rowname %in% insignificant_features, 0, 1))

  plot_data$rowname <- factor(plot_data$rowname, levels = plot_order)

  return(plot_data)
}

plot_dbrda <- function(plot_data) {
  # Plot the data as a horizontal bar plot
  g <- ggplot(data = plot_data, aes(x = value, y = rowname2, fill = variable)) +
    geom_bar(aes(alpha = significant), stat = "identity", position = position_dodge2()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    xlab("Effect size") +
    ylab("Feature") +
    labs(fill = "") + # Hide title of the legend
    guides(alpha = "none") + # Hide alpha legend
    scale_fill_manual(
      values = c("#60A68B", "#1F4068"),
      labels = c("Univariate", "Multivariate")
    ) +
    scale_alpha_continuous(range = c(0.5, 1), limits = c(0, 1))

  return(g)
}
