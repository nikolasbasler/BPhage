library(tidyverse)
library(patchwork)

alignment_plot <- function(blast_result) {
  blast_result %>%
    add_row(stitle="Query", sseqid="", qstart=1, qend = blast_result$qlength[1]) %>% # Add dummy row
    slice(c(n(), 1:(n() - 1))) %>% # Move dummy row to the bottom of the dataframe
    mutate(text = paste(stitle, sseqid)) %>%
    mutate(text = ifelse(nchar(text)>150, sub("^(.{0,})(.{145})$", "[...]\\2", text), text)) %>% # Shorten overly-long labels
    mutate(text = factor(text, levels = rev(unique(text)))) %>%
    ggplot(aes(x = qstart, xend = qend, y = text, yend = text, fill = pident)) +
    geom_segment(aes(color = pident), linewidth = 2) +
    labs(title = unique(blast_result$qseqid), x = "Position", y = NULL, color = "percent\nidentity\n") +
    scale_colour_gradient(limits = c(0, 100), na.value = "black")
}

# Load in data
families <- list.files("output/extra_searches") %>%
  sub("_contigs.*", "", .) %>%
  unique()

blastn <- list()
diamondx <- list()
for (f in families) {
  temp <- read.delim(paste0("output/extra_searches/",f,
                            "_contigs.blastn.tsv")) %>%
    mutate(qlength = as.numeric(str_split_i(qseqid, "_", 4)))
  blastn[[f]] <- split(temp, temp$qseqid,)
  
  temp <- read.delim(paste0("output/extra_searches/",f,
                            "_contigs.diamondx.tsv")) %>%
    mutate(qlength = as.numeric(str_split_i(qseqid, "_", 4)))
  diamondx[[f]] <- split(temp, temp$qseqid)
}

system("mkdir -p output/R/blast_and_diamond_plots")

# BLAST plots
for (f in families) {
  if (length(blastn[[f]])==0) {
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = paste0("No BLAST hits for\n",f))) +
      theme_void() 
    ggsave(paste0("output/R/blast_and_diamond_plots/", f, ".blast.NO.HITS.pdf"),
           empty_plot, width=5, height=0.5)
    next
  }
  plotlist <- list()
  hits <- c()
  longest_texts <- c()
  for (contig in names(blastn[[f]])) {
    plotlist[[contig]] <- alignment_plot(blastn[[f]][[contig]])
    hits <- blastn[[f]][[contig]]$sseqid %>%
      unique() %>%
      length() %>%
      c(hits, .) %>%
      subset(., .!=0)
    longest_texts <- paste(blastn[[f]][[contig]]$stitle, blastn[[f]][[contig]]$stitle) %>%
      nchar() %>%
      max() %>%
      c(longest_texts, .)
  }
  wrap <- wrap_plots(plotlist, ncol=1, heights=hits, guides="collect")
  total_hits <- sum(hits)
  ggsave(paste0("output/R/blast_and_diamond_plots/", f, ".blastn.pdf"), wrap, 
         width = 7.5 + min(150, max(longest_texts))/15, height = length(plotlist) + total_hits/5,
         limitsize = FALSE)
}

# DIAMOND plots
for (f in families) {
  if (length(diamondx[[f]])==0) {
    contig_names <- paste(names(diamondx[[f]]), collapse="\n")
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = paste0("No Diamond hits for\n",f))) +
      theme_void() 
    ggsave(paste0("output/R/blast_and_diamond_plots/", f, ".diamond.NO.HITS.pdf"),
           empty_plot, width=5, height= 0.5)
    next
  }
  plotlist <- list()
  hits <- c()
  longest_texts <- c()
  for (contig in names(diamondx[[f]])) {
    plotlist[[contig]] <- alignment_plot(diamondx[[f]][[contig]])
    hits <- diamondx[[f]][[contig]]$sseqid %>%
      unique() %>%
      length() %>%
      c(hits, .) %>%
      subset(., .!=0)
    longest_texts <- paste(diamondx[[f]][[contig]]$stitle, diamondx[[f]][[contig]]$stitle) %>%
      nchar() %>%
      max() %>%
      c(longest_texts, .)
  }
  wrap <- wrap_plots(plotlist, ncol=1, heights=hits, guides="collect")
  total_hits <- sum(hits)
  ggsave(paste0("output/R/blast_and_diamond_plots/", f, ".diamond.pdf"), wrap, 
         width = 7.5 + min(150, max(longest_texts))/15, height = length(plotlist) + total_hits/5,
         limitsize = FALSE)
  }
