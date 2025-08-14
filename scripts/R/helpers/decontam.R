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

  contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.05)

  
  print(paste(sum(contamdf.prev$contaminant),"contaminant contigs detected:"))
  contamdf.prev %>% 
    filter(contaminant) %>%
    rownames_to_column("contig") %>%
    left_join(., classification, by="contig") %>%
    print()
  
  contaminants <- contamdf.prev %>%
    rownames_to_column("contig") %>%
    filter(contaminant) %>%
    select(contig) %>%
    unlist(use.names = FALSE) 
  
  return(list(contaminants=contaminants, plot=p))
    
}