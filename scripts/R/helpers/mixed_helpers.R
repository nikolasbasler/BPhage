forest_plot <- function(tbl, axis_name = NULL, plot_title = NULL) {
  tbl %>%
    # ggplot(aes(x = Estimate, y = axis_labels)) +
    # geom_pointrange(aes(xmin = Estimate - `Std. Error`,
    #                     xmax = Estimate + `Std. Error`,
    ggplot(aes(x = estimate, y = axis_labels)) +
    geom_pointrange(aes(xmin = estimate - error,
                        xmax = estimate + error,
                        alpha = p_adjusted <= 0.05),
                    color = "black",
                    size = 0.7) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    geom_text_repel(aes(label = sprintf("%.3f", p_adjusted),
                        color = p_adjusted <= 0.05),
                    direction = "y",
                    nudge_y = 0.1,
                    segment.size = 0,
                    box.padding = 0.5) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(y = axis_name,
         x = "Estimated slope") +
    ggtitle(plot_title)
}


logistic_fun <- function(x, intercept, slope) {
  1 / (1 + exp(-(intercept + slope * x)))
}

logistic_effect_fun <- function(s, h, l) {
  odds_ratio <- exp(s * (h - l))
  return(odds_ratio)
}

linear_fun <- function(x, intercept, slope) {
  intercept + slope * x
}

linear_effect_fun <- function(s, h, l) {
  log_change <- s * (h - l)
  fold_change <- 10^log_change
  return(fold_change)
}

ct_effect_fun <- function(s, h, l) {
  ct_change <- s * (h - l)
  return(ct_change)
}

mixed_model_plot <- function(filt_test_tibble, transform_fun, effect_fun, dark_col, bright_col, y_axis_label) {
  
  x_axis_label <- str_replace_all(filt_test_tibble$Item, "_", " ")
  inter <- filt_test_tibble$intercept
  inter_sd <- filt_test_tibble$sd_intercept
  slo <- filt_test_tibble$Estimate
  slo_sd <- filt_test_tibble$`Std. Error`
  high <- filt_test_tibble$highest
  low <- filt_test_tibble$lowest

  y_of_high <- transform_fun(high, intercept = inter, slope = slo)
  y_of_low <- transform_fun(low, intercept = inter, slope = slo)
  
  effect_size <- effect_fun(s = slo, h = high, l = low)

  plot_min <- low - (high - low) / 5
  plot_max <- high + (high - low) / 5
  
  line_tibble <- tibble(x_val = seq(plot_min, plot_max, length.out = 100),
                        y_val = transform_fun(x_val, intercept = inter, slope = slo),
                        y_lower = transform_fun(x_val, intercept = inter - inter_sd, slope = slo - slo_sd),
                        y_upper = transform_fun(x_val, intercept = inter + inter_sd, slope = slo + slo_sd))
  
  y_stretched <- c(NA, NA)
  # if (slo < 0 ) {
  if (filt_test_tibble$which_y_end_to_stretch == "lower_end") {
    y_stretched <- c(min(line_tibble$y_lower) * filt_test_tibble$y_stretching_factor, NA)
  }
  if (filt_test_tibble$which_y_end_to_stretch == "upper_end") {
    y_stretched <- c(NA, max(line_tibble$y_upper) * filt_test_tibble$y_stretching_factor)
  }
  
  # Arrow Coordinates
  arrow_x_center <- (low + high) / 2
  arrow_width <- (high - low) / 15  # Adjust to control how "fat" the arrow is
  
  arrow_base <- y_of_low
  arrow_tip <- y_of_high
  
  # More coordinates (used later. Yes, it's a bit over-engineered...)
  arrow_head_base_y <- arrow_tip - (arrow_tip - arrow_base) * 0.2
  arrow_head_center_y <- (arrow_tip + arrow_head_base_y) / 2
  relative_arror_length <- abs(arrow_head_base_y - arrow_base) /
    max(line_tibble$y_upper - min(line_tibble$y_lower))
  
  # Define polygon for vertical arrow
  arrow_poly <- tibble(
    x = c(
      arrow_x_center - arrow_width,  # bottom left
      arrow_x_center + arrow_width,  # bottom right
      arrow_x_center + arrow_width,  # shaft right
      arrow_x_center + 2 * arrow_width,  # arrowhead right
      arrow_x_center,               # arrow tip
      arrow_x_center - 2 * arrow_width,  # arrowhead left
      arrow_x_center - arrow_width,  # shaft left
      arrow_x_center - arrow_width   # back to bottom left (close shape)
    ),
    y = c(
      arrow_base,  # bottom
      arrow_base,  # bottom
      arrow_tip - (arrow_tip - arrow_base) * 0.2,  # shaft top
      arrow_tip - (arrow_tip - arrow_base) * 0.2,  # arrowhead base right
      arrow_tip,  # tip
      arrow_tip - (arrow_tip - arrow_base) * 0.2,  # arrowhead base left
      arrow_tip - (arrow_tip - arrow_base) * 0.2,  # shaft top
      arrow_base   # close
    )
  )
  
  mmplot <- ggplot(line_tibble, aes(x = x_val, y = y_val)) +
    
    # Draw ribbon for standard error 
    geom_ribbon(aes(ymin = y_lower, ymax = y_upper), fill = bright_col, alpha = 0.6) +
    
    # Dashed lines from low point to arrow
    annotate(
      "segment",
      # x = low,
      x = -Inf,
      xend = arrow_x_center,
      y = arrow_base,
      yend = arrow_base,
      linetype = "dashed",
      color = "black"
    ) +
    annotate(
      "segment",
      x = high,
      # x = -Inf,
      # xend = arrow_x_center,
      xend = -Inf,
      y = arrow_tip,
      yend = arrow_tip,
      linetype = "dashed",
      color = "black"
    ) +
    
    # Draw curve and points for min and max observed values
    geom_line() +
    geom_point(x = high, y = y_of_high, size = 4) +
    geom_point(x = low, y = y_of_low, size = 4) +
    
    # Draw the arrow
    geom_polygon(data = arrow_poly, aes(x = x, y = y), fill = dark_col, alpha = 1) +
    
    # Print the effect size into the arrow
    annotate(
      "text",
      x = arrow_x_center,
      y = (arrow_base + arrow_head_base_y) / 2,
      # label = paste0("OR = ", round(effect_size, 2)),
      label = sprintf("%.2f", effect_size),
      angle = 90,
      # size = ifelse(relative_arror_length < 0.35, 2.75, 4),
      size = 4,
      color = "white"
    ) +
    
    # Print asterisks for siginificance level into head of arrow
    annotate(
      "text",
      x = arrow_x_center,
      y = arrow_head_center_y,
      label = filt_test_tibble$p_adjust_significant,
      color = "white",
      vjust = ifelse(slo < 0, 0.5, 1),
      # size = ifelse(relative_arror_length < 0.35, 3.5, 5),
      size = 4.5,
      # fontface = "bold"
    ) +
    
    # Make axis ticks only for the coordinates of the points
    scale_x_continuous(
      breaks = c(low, high),
      labels = c(round(low, 0), round(high, 0))) +
    scale_y_continuous(
      limits = y_stretched, # This will stretch the y-axis, i.e. compress the arrow
      breaks = c(y_of_low, y_of_high),
      labels = c(round(y_of_low, 2), round(y_of_high, 2))) +
    
    # Axis texts, title and theme
    labs(x = x_axis_label, y = y_axis_label) +
    # ggtitle(filt_test_tibble$gene) +
    
    
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
      ) 
    
  return(mmplot)
  
}


legend_factory <- function(title, items, colors, position) {
  legend_df <- tibble(items = factor(items, levels = items),
                      x = seq_along(items),
                      y = 1)
  
  dummy_plot <- ggplot(legend_df, aes(x, y, color = items)) +
    # draw zero‐size points (so nothing appears in the panel)
    geom_point(size = -1, show.legend = TRUE) +
    
    # define exactly the keys you want
    scale_color_manual(
      name   = title,
      values = setNames(colors, items),
      breaks = items
    ) +
    
    # force the legend symbols to show at a normal size
    guides(color = guide_legend(
      override.aes = list(size = 5)
    )) +
    
    # strip out axes / grid
    theme_void() +
    theme(
      legend.position = position
    )
  
  legend_only <- dummy_plot + 
    guide_area() +
    plot_layout(
      # define a 2‑row layout: "A" is the plot, "L" is the legend
      design  = "A\nL",
      guides  = "collect",
      # give the plot row zero height so only the legend shows
      heights = c(0, 1)
    )
  return(legend_only)
}



# 
# layered_p_adjustments <- function(slop = slopes, gene_or_pathogen = "gene") {
#   
#   p_value_list <- list()
#   hypotheses <- c()
#   for (level in names(slop)) {
#     p_value_list[[level]] <- slop[[level]] %>%
#       arrange(test_name) %>%
#       select(test_name, raw_p_value) %>%
#       deframe()
#     hypotheses <- c(hypotheses, names(p_value_list[[level]]))
#   }
#   
#   layers <- list()
#   layers$genes <- slop$cropland[[gene_or_pathogen]]
#   layers$cropland <- "Cropland_in_2km_radius"
#   layers$pesticides <- "Pesticides (total)"
#   layers$pest_groups <- unique(slop$pest_groups$Item)
#   layers$specific$`Fungicides and Bactericides` <- slop$specific_pests %>%
#     filter(str_detect(Item, "Fung & Bact ")) %>%
#     distinct(Item) %>%
#     unlist(use.names = FALSE)
#   layers$specific$Herbicides <- slop$specific_pests %>%
#     filter(str_detect(Item, "Herbicides ")) %>%
#     distinct(Item) %>%
#     unlist(use.names = FALSE)
#   layers$specific$Insecticides <- slop$specific_pests %>%
#     filter(str_detect(Item, "Insecticides ")) %>%
#     distinct(Item) %>%
#     unlist(use.names = FALSE)
#   
#   # Define the initial weights.
#   # Here we allocate alpha equally to the Stage 1 tests (genes) and zero to the rest.
#   # When the hypothesis of a gene-test is rejected, its weight will be transferred.
#   w <- rep(0, length(hypotheses))
#   names(w) <- hypotheses
#   w[names(p_value_list$cropland)] <- 1 / length(p_value_list$cropland)
#   
#   # Create an empty transition matrix.
#   transitions <- matrix(0, nrow = length(hypotheses), ncol = length(hypotheses),
#                         dimnames = list(hypotheses, hypotheses))
#   for (gene in layers$genes) {
#     # Layer 1 to 2: cropland fraction to total pest
#     L1_test_name <- paste0(gene, "; Cropland_in_2km_radius")
#     L2_test_name <- paste0(gene, "; Pesticides (total)")
#     
#     if ( L1_test_name %in% rownames(transitions) & L2_test_name %in% colnames(transitions) ) {
#       transitions[L1_test_name, L2_test_name] <- 1
#     }
#     
#     # Later 2 to 3: total pest to pest group
#     for (Layer3 in layers$pest_groups) {
#       L3_test_name <- paste0(gene, "; ", Layer3)
#       
#       if ( L2_test_name %in% rownames(transitions) & L3_test_name %in% colnames(transitions) ) {
#         transitions[L2_test_name, L3_test_name] <- 1 / length(layers$pest_groups)
#       }
#       
#       # Layer 3 to 4: pest group to specitic pest
#       for (Layer4 in layers$specific[[Layer3]]) {
#         L4_test_name <- paste0(gene, "; ", Layer4)
# 
#         if ( L3_test_name %in% rownames(transitions) & L4_test_name %in% colnames(transitions) ) {
#           transitions[L3_test_name, L4_test_name] <- 1 / length(layers$specific[[Layer3]])
#         }
#       }
#     }
#   }
#   # Read the transition matrix line by line, where the row name is the test in 
#   # question and the column name is the test in the next layer of testing.
#   # A 0 means that this test in the row name is not propagted to the test in the 
#   # col name. A number different to 0 is the weight this test (row name) will 
#   # have on the test in the next layer (col name). The weights in each line 
#   # should sum up to 1. In an easy setting, equally divide the weight to the
#   # test of the next layer.
#   
#   graph <- matrix2graph(m = transitions, weights = w)
#   
#   p_raw <- c(p_value_list$cropland, p_value_list$total_pest, p_value_list$pest_groups, p_value_list$specific_pests)
#   # Optionally, inspect the graph.
#   # print(graph)
#   
#   # Cor mat
#   cor_mat <- transitions %>%
#     as.data.frame() %>%
#     rownames_to_column("test_name") %>%
#     pivot_longer(-test_name) %>%
#     mutate(value = ifelse(value != 0, NA, 0),
#            value = ifelse(test_name == name, 1, value)) %>%
#     pivot_wider() %>%
#     column_to_rownames("test_name") %>%
#     as.matrix()
#   # Ensure the matrix is symmetric by copying the upper triangle to the lower triangle
#   # (if your procedure requires a symmetric correlation matrix)
#   cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)]
# 
#   result <- gMCP(pvalues = p_raw, graph = graph)
# 
#   adjusted_p <- tibble(test_name = names(result@adjPValues),
#                        p_adjusted = result@adjPValues)
#   
#   subgraphs <- transitions %>%
#     as.data.frame() %>%
#     rownames_to_column("test_name") %>%
#     as_tibble() %>%
#     pivot_longer(-test_name, names_to = "child") %>%
#     filter(value > 0) %>%
#     select(-value) %>%
#     mutate(gene = str_split_i(test_name, ";", 1),
#            parent_layer = str_split_i(test_name, "; ", 2),
#            child_layer = str_split_i(child, "; ", 2),
#     ) %>%
#     select(gene, parent_layer, child_layer) %>%
#     distinct()
#   
#   graph_plots <- list()
#   for (gen in unique(subgraphs$gene)) {
#     for (parent in unique(subgraphs$parent_layer)) {
#       children <- subgraphs %>%
#         filter(gene == gen,
#                parent_layer == parent) %>%
#         select(child_layer) %>%
#         unlist(use.names = FALSE)
#       
#       filt_names <- tibble(test_name = colnames(transitions)) %>%
#         separate_wider_delim(test_name, delim = "; ", names = c("gene", "test")) %>%
#         filter(gene == gen,
#                test %in% parent | test %in% children) %>%
#         mutate(test_name = paste0(gene, "; ", test)) %>%
#         select(test_name) %>%
#         unlist(use.names = FALSE)
#       
#       filtered_matrix <- transitions %>%
#         as.data.frame() %>%
#         rownames_to_column("test_name") %>%
#         tibble() %>%
#         filter(test_name %in% filt_names) %>%
#         select(test_name, all_of(filt_names)) %>%
#         column_to_rownames("test_name") %>%
#         as.matrix()
#       
#       if (nrow(filtered_matrix) > 0 ) {
#         graph_plots[[parent]][[gen]] <- hGraph(m = filtered_matrix, nHypotheses = nrow(filtered_matrix), nameHypotheses = rownames(filtered_matrix))
#       }
#     }
#   }
#   
#   return(list(transition_matrix = transitions, gMCP_result = result, subgraph_plots = graph_plots, adjusted_p_values = adjusted_p))
# }
