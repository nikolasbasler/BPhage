# BPhage data analysis
## Versions
- Viper version 2.1 (commit 8244c2eeeebeb2c1fbe4082034ca329e6f0a4a10 13 March 2024, https://github.com/Matthijnssenslab/ViPER)
- For software versions see [data/env_viper_bphage.yml](data/env_viper_bphage.yml) and [data/env_ncbi.yml](data/env_ncbi.yml)
- R 4.3.1
- RStudio 2023.12.1+402 (2023.12.1+402)
- To reproduce R package versions run `renv::restore()`

## Script order
### Assembly
- Prepare: Copy raw files from K-drive. Sample 4295 (BE_16561_aut_rec_d) failed in sequencing run R4317 and was re-run in R4341.
- `scripts/HPC/rename.BPhage.nucleomics.sh` (contained in R4317 folder on K-drive)
    - Requires: `rename.BPhage.nucleomics.scheme` (also contained in R4317)
    - Output: Properly named raw fastq files (not by tube ID, as before, but by sample ID).
- Download bee ref seq and index with bowtie2: `scripts/HPC/download_bee_ref_genome.slrm`
- ViPER: `scripts/HPC/bphage_viper_with_dedup.slrm` (array of 471)
    - Requires: 
        - `data/BPhage.sample.list`
        - Raw read data in `$VSC_SCRATCH/BPhage/raw`
    - Output: Trimmed reads, trimmed host-removed reads, assembly. No Krona (done before with previous version).
- Re-organise ViPER output: `scripts/HPC/reorganise_viper_output.slrm`
    - Requires: 
        - `data/BPhage.sample.list`
        - ViPER output
    - Output: Oririnal files from each sample's viper output are moved into a common `CONTIGS`, `QC/FASTQC` (also with multiqc) `QC/QUAST` and `READ` (containing deduped trimmed and hostout reads) folder inside `output/bphage_viper_output`
- Cross sample clustering: `scripts/HPC/cross_sample_clustering.slrm`
    - Requires: Re-organised ViPER assemblies in `output/bphage_viper_output/CONTIGS/`
    - Output: Cross-sampled file: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`, `output/bphage_ALL_1kb_cross_95-85_clusters.tsv.gz`

### Phage identification
- CheckV (can be done on parallel with geNomad and mapping): `scripts/HPC/checkv.slrm`
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_checkv`
- geNomad (can be done in parallel with CheckV and mapping): `scripts/HPC/genomad.slrm`
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_genomad`
- Filter for >= 70% complete phages: `scripts/HPC/filter_classification.slrm`
    - Requires: 
        - Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
        - CheckV output
        - geNomad output
        - ICTV's VMR table: `data/VMR_19-250422_MSL37.xlsx` (downloaded from https://ictv.global/vmr on 14 March 2024)
        - Python script to do the actual filtering: `scripts/HPC/filter_classification.py` (called by the slrum script)
    - Output: 
        - Fasta file and merged genomad and checkv table with >= 50%-complete phages: `output/bphage_ALL_1kb_phages`
        - Fasta file and merged genomad and checkv table with >= 50%-complete unclassified viruses (containing contigs that are either "Unclassified" or classified only as "Viruses"): `output/bphage_ALL_1kb_unclassified_viruses`

### Mapping
- Prepare mapping: `scripts/HPC/prepare_mapping.slrm`
    - Requires: 
        - Phage contigs at `output/bphage_ALL_1kb_phages.fasta.gz`
        - Unclassified virus contigs at `output/bphage_ALL_1kb_unclassified_viruses.fasta.gz`
        - Genomad summary to get a list of Picobirnaviridae contigs: `output/bphage_ALL_1kb_genomad/bphage_ALL_1kb_cross_95-85_summary/bphage_ALL_1kb_cross_95-85_virus_summary.tsv.gz`
        - Entire assembly to add Picornaviridae contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: bwa-indexed mapping ref containing all contigs identified as phage, unclassified viruses and Picovridae contigs: `$VSC_SCRATCH/BPhage/ref/bphage_mapping_ref.fasta`
- Mapping: `scripts/HPC/bphage_mapping.slrm` (array of 471)
    - Requires: 
        - `data/BPhage.sample.list`
        - bwa-indexed mapping ref at `$VSC_SCRATCH/BPhage/ref/bphage_mapping_ref.fasta`
    - Output: 
        - Mapping alignments in `output/mapped`
        - Mapping stats in `output/mapped`
- Gather mapping stats: `scripts/HPC/gather_mapping_stats.slrm`
    - Requires:
        - `data/BPhage.sample.list`
        - Mapping stats of individual samples in `output/mapped/`
        - Mapping ref to get the complete contig list: `$VSC_SCRATCH/BPhage/ref/bphage_mapping_ref.fasta`
    - Output:
        - Separate tables for mapped reads, horizontal coverage and mean depth (samples in columns, contigs in rows) and a table with all contig lengths in `output/mapping_stats`
    