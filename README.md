# BPhage data analysis
## Versions
- Viper
- etc.

## 
## Script order
- Prepare: Copy raw files from K-drive. Sample 4295 (BE_16561_aut_rec_d) failed in sequencing run R4317 and was re-run in R4341.
- `rename.BPhage.nucleomics.sh` (contained in R4317 folder on K-drive)
    - Requires `rename.BPhage.nucleomics.scheme` (also contained in R4317)
    - Output: Fastq files were named based on tube ID, renmane them based on sample ID.
- ViPER: `bphage_viper_with_dedup.slrm` (array of 471, so run in 4 batches: 3x118 1x117)
    - Requires: `data/BPhage.sample.list`, nr database `nr_20231025`