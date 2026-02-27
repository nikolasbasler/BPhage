forest_plot <- function(tbl, axis_name = NULL, plot_title = NULL) {
  tbl %>%
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

extract_legend <- function(ggplot_object) {
  extracted_legend <- ggplot_object  + 
    guide_area() +
    plot_layout(
      guides  = "collect",
      heights = c(0, 1))
  return(extracted_legend)
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
  inter_sd <- filt_test_tibble$se_intercept
  high <- filt_test_tibble$highest
  low <- filt_test_tibble$lowest
  slo <- filt_test_tibble$Estimate
  y_of_high <- transform_fun(high, intercept = inter, slope = slo)
  y_of_low <- transform_fun(low, intercept = inter, slope = slo)
  slo_sd <- ifelse(high < 0 , -filt_test_tibble$`Std. Error`, filt_test_tibble$`Std. Error`)

  effect_size <- effect_fun(s = slo, h = high, l = low)
  effect_prefix <- ""
  effect_suffix <- ""
  
  # For percent change rather than fold change in gene relabund plots:
  if (y_axis_label == "Log rel. gene abund.") {
    effect_size <- round((effect_size - 1) * 100)
    effect_suffix <- "%"
    if ( effect_size >= 1) {
      effect_prefix <- "+"
    }
  } else {
    effect_size <- sprintf("%.2f", effect_size)
  }
  
  plot_min <- low - (high - low) / 5
  plot_max <- high + (high - low) / 5
  
  line_tibble <- tibble(x_val = seq(plot_min, plot_max, length.out = 100),
                        y_val = transform_fun(x_val, intercept = inter, slope = slo),
                        y_lower = transform_fun(x_val, intercept = inter - inter_sd, slope = slo - slo_sd),
                        y_upper = transform_fun(x_val, intercept = inter + inter_sd, slope = slo + slo_sd))

  y_stretched <- c(NA, NA)
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
    # annotate(
    #   "text",
    #   x = arrow_x_center,
    #   y = (arrow_base + arrow_head_base_y) / 2,
    #   # label = sprintf("%.2f", effect_size),
    #   label = paste0(effect_prefix, effect_size, effect_suffix),
    #   angle = 90,
    #   # size = 4,
    #   size = 3, 
    #   color = "white"
    # ) +
    
    # Print asterisks for siginificance level into head of arrow
    # annotate(
    #   "text",
    #   x = arrow_x_center,
    #   y = arrow_head_center_y,
    #   label = filt_test_tibble$p_adjust_significant,
    #   color = "white",
    #   # vjust = ifelse(slo < 0, 0.5, 1),
    #   vjust = ifelse(slo < 0, 0.75, 1),
    #   size = 4.5,
    #   # fontface = "bold"
    # ) +
    
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
      legend.position = position,
      legend.title = element_text(face = "bold")
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

diagnostics_logistic_model <- function(model, name = "") {
  simRes <- simulateResiduals(fittedModel = model, n = 1000)
  resDat <- data.frame(
    predicted = simRes$fittedPredictedResponse,
    residual  = simRes$scaledResiduals
  )
  
  p <- ggplot(resDat, aes(x = predicted, y = residual)) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(
      x    = "Predicted probability",
      y    = "DHARMa scaled residuals",
      title= name
    ) +
    theme_minimal()
  return(p)
}

diagnostics_linear_model <- function(model, name = "") {
  resDat <- data.frame(
    fitted   = fitted(model),
    residual = resid(model)
  )
  
  p <- ggplot(resDat, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x    = "Fitted values",
      y    = "Raw residuals",
      title= name
    ) +
    theme_minimal()
  return(p)
}