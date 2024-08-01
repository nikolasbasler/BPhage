
##  LOAD lines 1-77 from alphabeta script

bacphlip <- read.delim("output/lifestyle/bacphlip/bphage.fasta.bacphlip") %>%
  rename("contig" = "X")
View(bacphlip)
