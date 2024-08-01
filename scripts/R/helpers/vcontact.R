

add_dataset_to_cytoscape_table <- {
  ## ADD COLUMN TO DATA
  ## - Export node list from cytoscape
  ## - Open it here
  ## - Add column containing the dataset
  ## - Save it and import into cytoscape
  
  
  final_assignments <- read.csv("output/vcontact3/final_assignments.csv") %>%
    select(-X) %>%
    select(RefSeqID, Proteins, Reference, contains("Size"), contains("prediction"), network)
  
  file_location <- "~/Library/CloudStorage/OneDrive-KULeuven/PhD/Virome/BPhage/output/vcontact3/"
  bins <- c("all", "bin_0", "bin_1", "bin_2", "bin_3")
  
  cyjs <- list()
  result_table <- list()
  for (bin in bins) {
    if (bin == "all") {
      cyjs[[bin]] <- read.csv(paste0(file_location, "graph.cyjs default node.csv"))
    } else {
      cyjs[[bin]] <- read.csv(paste0(file_location, "graph.", bin,".cyjs default node.csv"))
    }
    cyjs[[bin]] %>% 
      mutate_all(~ifelse(is.na(.), "", .)) %>%
      mutate(dataset = "vcontact3") %>%
      mutate(dataset = ifelse(str_detect(name, "NODE_"), "BPhage", dataset)) %>%
      mutate(dataset = ifelse(str_detect(name, "Bonilla"), "BonillaRosso", dataset)) %>%
      mutate(dataset = ifelse(str_detect(name, "Deboutte"), "Deboutte", dataset)) %>%
      mutate(dataset = ifelse(str_detect(name, "Busby"), "Busby", dataset)) %>%
      left_join(., final_assignments[], by = join_by(name == RefSeqID)) %>%
      write_csv(paste0(file_location, "graph.", bin, ".cyjs.with_dataset_column.csv"), quote = "all")
  } 
}


