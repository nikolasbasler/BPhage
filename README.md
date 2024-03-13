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
- Download bee ref seq and index with bowtie2: `download_bee_ref_genome.slrm`
- ViPER: `bphage_viper_with_dedup.slrm` (array of 471, so run in batches)
    - Requires: 
        - `data/BPhage.sample.list`
        - Raw read data in `$VSC_SCRATCH/BPhage/raw`
    - Output: Trimmed reads, trimmed host-removed reads, assembly. No Krona (done before with previous version).
- Re-organise ViPER output: `reorganise_viper_output.slrm`
    - Requires: `data/BPhage.sample.list`
    - Output: Oririnal files from each sample's viper output are moved into a common `CONTIG`, `QC/FASTQC` (also with multiqc) `QC/QUAST` and `READ` (containing deduped trimmed and hostout reads) folders. 