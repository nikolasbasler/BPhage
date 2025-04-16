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
#   layers$cropland <- "ha_cropland_in_2k_radius"
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
#     L1_test_name <- paste0(gene, "; ha_cropland_in_2k_radius")
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
