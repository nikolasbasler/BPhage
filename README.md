# BPhage data analysis
## Versions
- Viper version 2.1 (commit 8244c2eeeebeb2c1fbe4082034ca329e6f0a4a10 13 March 2024)
- For software versions see [data/viper_bphage_env.yml](data/viper_bphage_env.yml)

## 
## Script order
- Prepare: Copy raw files from K-drive. Sample 4295 (BE_16561_aut_rec_d) failed in sequencing run R4317 and was re-run in R4341.
- `rename.BPhage.nucleomics.sh` (contained in R4317 folder on K-drive)
    - Requires: `rename.BPhage.nucleomics.scheme` (also contained in R4317)
    - Output: Properly named raw fastq files (not by tube ID, as before, but by sample ID).
- Download bee ref seq and index with bowtie2: `download_bee_ref_genome.slrm`
- ViPER: `bphage_viper_with_dedup.slrm` (array of 471, so run in batches)
    - Requires: 
        - `data/BPhage.sample.list`
        - Raw read data in `$VSC_SCRATCH/BPhage/raw`
    - Output: Trimmed reads, trimmed host-removed reads, assembly. No Krona (done before with previous version).
- Re-organise ViPER output: `reorganise_viper_output.slrm`
    - Requires: 
        - `data/BPhage.sample.list`
        - ViPER output
    - Output: Oririnal files from each sample's viper output are moved into a common `CONTIGS`, `QC/FASTQC` (also with multiqc) `QC/QUAST` and `READ` (containing deduped trimmed and hostout reads) folder inside `output/bphage_viper_output`
- Cross sample clustering: `cross_sample_clustering.slrm`
    - Requires: Re-organised ViPER assemblies in `output/bphage_viper_output/CONTIGS/`
    - Output: Cross-sampled file: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`, `output/bphage_ALL_1kb_cross_95-85_clusters.tsv.gz`
- CheckV (can be done on parallel with geNomad and mapping): `checkv.slrm`
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_checkv`
- geNomad (can be done in parallel with CheckV and mapping): `genomad.slrm`
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_genomad`
- Filter for >= 50% complete phages: `filter_classification.slrm`
    - Requires: 
        - Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
        - CheckV output
        - geNomad output
        - ICTV's VMR table: `data/VMR_MSL38_v2.csv` (derived from `VMR_MSL38_v2.xlsx`, downloaded from https://ictv.global/vmr on 14 March 2024)
        - R script to do the actual filtering: `filter_classification.R` (called by the slrum script)
    - Output: 
        - Fasta file and merged genomad and checkv table with >= 50%-complete phages
        - Fasta file and merged genomad and checkv table with >= 50%-complete unclassified viruses (containing contigs that are either "Unclassified" or classified only as "Viruses")
- Mapping (can be done in parallel with geNomad and CheckV):