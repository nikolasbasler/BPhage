decontam_identification <- function(abtable, classification) {
  # Returns a vector with contaminants. And a plot Only the "prevalence" method 
  # is used, with a more stringent p-value threshold of 0.01.
  # The removed contigs will be reported for review. 
  # The script will return the abundance table without the contaminants and 
  # without the blank controls.
  
  # Removing samples with only zeros. Will also be done in main script.
  abtable_filt <- abtable[rowSums(abtable) > 0, ]
  abtable_filt <- abtable[, colSums(abtable) > 0]
  
  # Make Phyloseq object
  otu = otu_table(abtable_filt, taxa_are_rows = TRUE)
  ps_metadata = sample_data(metadata)
  ps = phyloseq(otu, ps_metadata)
 
  #### Creating a dataframe for plotting
  df = as.data.frame(sample_data(ps))
  df$LibrarySize = sample_sums(ps) # This sums up the total number of mapped reads for each sample.
  df <- df[order(df$LibrarySize),] # This is giving a rank number based on the total number of mapped reads.
  df$Index <- seq(nrow(df)) # This is effectively ordering the samples by their total number of mapped reads.
  
  p <- df %>%
    ggplot(aes(x=Index, y=LibrarySize, color=Sample_or_control, label=Sample_ID)) + 
    geom_point()

  # Blank pool 4 has to be removed for this analysis, because it would otherwise
  # mess up the detection of contaminants.
  ps = subset_samples(ps, Sample_ID!="Blank_pool_04" )
  
  ############################################################
  #### Contamination identification by prevalence
  ############################################################
  
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "sample_blank"
  # sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "extraction_blank"
  # sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "extraction_blank" | sample_data(ps)$Sample_or_control == "sample_blank"
  
  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.05)
  # contamdf.freq <- isContaminant(ps, method="frequency", neg="is.neg", threshold=0.05, conc="Average_Qubit")

  
  print(paste(sum(contamdf.prev$contaminant),"contaminant contigs detected:"))
  contamdf.prev %>% 
    filter(contaminant) %>%
    merge(., classification, by=0) %>%
    column_to_rownames("Row.names") %>%
    print()
  
  contaminants <- contamdf.prev %>%
    rownames_to_column("contig") %>%
    filter(contaminant) %>%
    select(contig) %>%
    unlist(use.names = FALSE) 
  
  return(list(contaminants=contaminants, plot=p))
    
  # colnames(abtable)[abtable["NODE_A3568_length_1426_cov_8.824314_Blank_pool_02",] > 0] # Only occurs in 2 blanks and no sample -> contaminant
  # colnames(abtable)[abtable["NODE_A1569_length_1169_cov_7.784799_FR_19773_spr_ile_d",] > 0] # Only occurs in 1 blank and 1 sample (not processed on the same day) -> contaminant
  # colnames(abtable)[abtable["NODE_A516_length_2179_cov_6.701713_FR_19773_spr_ile_d",] > 0] # Only occurs in 1 blank and 1 sample (not processed on the same day) -> contaminant
}

crosscontam_plot <- function(tpm_table, classification, tax) {
  
  original_contigs <- classification %>%
    filter(Species==tax) %>%
    select(contig) %>%
    unlist(use.names=FALSE)
  tpm_tax_blanks <- tpm_table %>%
    filter(Species==tax) %>%
    select(Species, contains("Blank"))
  print(paste0(tax," TPM in blanks:"))
  print(tpm_tax_blanks)
  
  colnames(tpm_table)[1]="name"
  col_levels <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12")
  
  pool_info <- read.csv("data/blank_pool_info.csv") %>%
    select(Sample_ID, Extraction_plate_row, Extraction_plate_col) %>%
    distinct()
  
  ## BLANKS MAY HAVE HIGH TPM BECAUSE THEY DONT HAVE ANYTHING ELSE...
  
  plate_plt <- tpm_table %>%
    filter(name==tax) %>%
    column_to_rownames("name") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    inner_join(., metadata, by="Sample_ID") %>%
    left_join(., pool_info, by=c("Sample_ID")) %>% 
    filter(Sample_or_control!="WTA_blank") %>%
    mutate(Extraction_plate_row = coalesce(Extraction_plate_row.x, Extraction_plate_row.y),
           Extraction_plate_col = coalesce(Extraction_plate_col.x, Extraction_plate_col.y)) %>%
    select(-c(Extraction_plate_row.x, Extraction_plate_col.y)) %>%
    select(Sample_ID, Country, Hive_ID, Season, Gut_part, Health, all_of(tax), Extraction_plate, Extraction_plate_row, Extraction_plate_col) %>%
    mutate(lab=paste0(paste(Country, Hive_ID, sep="_"),"\n", paste(Season, Gut_part, sep="_"),"\n",Health)) %>%
    mutate(lab = ifelse(grepl("NA", lab), Sample_ID, lab)) %>%
    mutate(Extraction_plate_row = factor(Extraction_plate_row, levels=rev(LETTERS[1:8]))) %>%
    mutate(Extraction_plate_col = factor(Extraction_plate_col, levels=col_levels)) %>%
    mutate(is_blank = grepl("Blank", Sample_ID)) %>%
    ggplot(aes(x=Extraction_plate_col, y=Extraction_plate_row, fill=.data[[tax]])) +
    geom_tile() +
    geom_text(aes(label = lab, color = is_blank), vjust = 0.5, size = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
    # geom_text(aes(label = lab), vjust = 0.5, color = "white", size=1) +
    facet_wrap(~Extraction_plate) +
    labs(x=NULL, y=NULL)
  
  return(plate_plt)
  
}
