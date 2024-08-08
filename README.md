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
- Download bee ref seq and and other datasets and index with bowtie2: `scripts/HPC/download_bee_ref_genome_and_additional_datasets.slrm`
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
    - Requires: 
        - Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
        - Data from other studies: `$VSC_SCRATCH/additional_datasets/other_studies.fa.gz`
    - Output: 
        - `output/bphage_ALL_1kb_checkv`
        - `output/other_studies_checkv`
- geNomad (can be done in parallel with CheckV and mapping): `scripts/HPC/genomad.slrm`
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_genomad`
- Filter for >= 50% complete phages: `scripts/HPC/filter_classification.slrm`
    - Requires: 
        - Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
        - CheckV output
        - geNomad output
        - ICTV's VMR table: `data/VMR_19-250422_MSL37.xlsx` (downloaded from https://ictv.global/vmr on 14 March 2024)
        - Python script to do the actual filtering: `scripts/HPC/filter_classification.py` (called by the slrum script)
    - Output: 
        - Fasta file and merged genomad and checkv table with >= 50%-complete phages, >= 50%-complete unclassified viruses (containing contigs that are either "Unclassified" or classified only as "Viruses") and Picobirnaviridae (no completeness threshold) for own dataset and for the additional datasets (Deboutte, Bonilla, Busby):
            - `output/bphage_ALL_1kb_phages.*`, `output/bphage_ALL_1kb_unclassified_viruses.*`, `output/bphage_picobirna.*`
            - `output/other_studies_phages.*`, `output/other_studies_unclassified_viruses.*`, `other_studies_picobirna.*` (the latter is empty because there are not Picobirnas in the additional datasets)
- Clustering with additional datasets: `clustering_with_additional_datasets.slrm`
    - Requires: Filtered assembly and contigs from additional datasets in `output/`
    - Output: File with cluster information, i.e. which contigs are 95%/85cov similar to bphage contigs: `output/bphage_and_others_clusters.tsv`

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

### Clustering
- vConTACT3: `scripts/HPC/vcontact3_clustering_with_inphared.slrm`
    - Requires: 
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
        - Phage and unclassified contigs from other studies: `output/other_studies_phages.fasta.gz`, `output/other_studies_unclassified_viruses.fasta.gz`
        - INPHARED dataset: `$VSC_SCRATCH/BPhage/additional_datasets/inphared_3Apr2024_genomes_excluding_refseq.fa.gz`
    - Output: `$repo_location/output/vcontact3/bphage_vcontact3_b38_with_inphared`

### Annotation
- Pharokka: `annotation_pharokka_bphage_and_others.slrm`
    - Requires:
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
        - Phage and unclassified contigs from other studies: `output/other_studies_phages.fasta.gz`, `output/other_studies_unclassified_viruses.fasta.gz` (no picobirnas found in the other studies)
    - Output: `output/output/annotation/pharokka/bpgage_and_others`
- Collabfold: Done by collaborator George Bouras.
    - Requires: Pharokka's prodigal-gv output: `output/annotation/pharokka/bpgage_and_others/prodigal-gv.faa` (in communication with George the file was named `bpgage_and_others_prodigal-gv.faa`)
    - Output: Protein structures (Alpha fold 2) of all identified proteins from bphage and other studies. Provided by George, placed into `output/annotation/phold_colabfold_structures`
- Prepare phold compare:
    - Requires: 
        - pdb files of George's collabfold structures in `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - List of contigs that were sent to George but are not in the bphage_and_others dataset anymore because they were either filtered out by geNomad/CheckV (from additional studies) or collapsed by clustering: `$repo_location/data/phold/sent_to_george_but_removed_from_original_dataset`
    - Output: 
        - Renamed pdb files with shortened names, otherwise phold will crash.
        - pdb files without corresponding entry in pharokka's output will be moved to `output/annotation/phold_colabfold_structures/basler_output_renamed/filtered_out_renamed_pdbs`
- Phold: `annotation_phold_compare.slrm`
    - Requires: 
        - Pharokka's gbk output: `output/annotation/pharokka_bphage_and_others/bphage_and_others.gbk`
        - Predicted structures at `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - sed script to re-lengthen contig names: `data/phold/relengthen.contig.names.sed`
    - Output:
        - Predicted functions at: `output/annotation/phold_compare_bphage_and_others`
