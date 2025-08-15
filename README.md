# The honey bee triad: Phages are mutualistic partners in the gut microbiome of *Apis mellifera*

This repository contains all scripts and usage instructions to reproduce the analysis for the paper by Basler et al (XXX ref). All 97 core phage genomes as well as all >90% complete non-core genomes (1047 additional sequences) are on Genbank (See Supplementary file, sheet i - Genbank accessions). The remaining 1199 non-core, <90% complete genomes can be found in this Github repository. 

This pipeline is split into two parts. The first part is meant for a high-performance computer (HPC) and can be skipped, if so wanted. The second part is for the statistical analysis and visualisation using RStudio. The scripts pretend to be on the same computer but it is possible to clone this repo to an HPC and a local computer, run the HPC scripts, copy the relevant output files from the HPC to the local computer and continue with the statistical analysis using the R project.

To clone the repository, please run
```
git clone https://github.com/nikolasbasler/BPhage
```
**Note**: The output of the tools and scripts will end up in the `output` folder inside the repo (which is why it's not tracked by git). The HPC scripts will create around 1.5 TB in total (plus intermediate storage, see below), the R scripts around 1 GB. Make sure to have enough free space or manage the output as it comes.

**If you want to skip the HPC part and only want to re-run the statistical analysis:** Clone the repo as mentioned above, extract the mid-save archive, and then skip ahead to the "R scripts" section of this README file.
```
tar -kxvzf mid_save.tar.gz
```

## HPC scripts
### General info
- All scripts for this section are located in `scripts/HPC`.
- These scripts mostly represent jobs for a slurm scheduler. 
- If you also use slurm, you will have to adapt the instructions at the beginning of each script according to your setup. Particularly the `account` name and probably also the resource allocation will be different on your system.
- If you want to port these instructions to a different scheduler, your favourite generative AI tool will be very helpful.
- You can also run all scripts without a scheduler, directly with `bash` but note that some of them run for a long time, especially if the computational resources are limited.
- Array scripts are iterations of the same slurm job, meant to run simultaneously. If you have to run them in sequence, you can embed the entire script into a loop, with `$SLURM_ARRAY_TASK_ID` as iterator, but this will take a **very** long time for most of them.

### Installations
For tool versions see the conda environment `.yml` files: `data/env_*.yml`
- Mamba/conda: Please follow the official documentation (e.g.: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). For this analysis, `mamba` v.1.4.2 and `conda` v.23.1.0 were used.
- ViPER pipeline: Please follow the instructions on Github **BUT** use `data/env_viper_bphage.yml` from this repository instead of ViPER's `viper.yml`: https://github.com/Matthijnssenslab/ViPER. For this analysis, ViPER version 2.1 was used.
- Additional conda environments: To install additional conda environments with all necessary tools, use the `env_*.yml` files in `data/`: 
- `env_cobra.yml`
- `env_iphop.yml`
- `env_mopup.yml`
- `env_ncbi.yml`
- `env_pharokka.yml`
- `env_phold.yml`
- `env_replidec.yml`
- `env_ska.yml`
- `env_vcontact3.yml`

E.g. like this:
```
mamba env create -f env_cobra.yml
```

- **Note**: The tool MOP-UP (used for taxonomic clustering of the microviruses) relies on a library called `Boost`, which I didn't manage to install properly via conda but instead relied on a pre-installed module. If it gives you trouble, see version information in `env_mop-up_boost_module.txt` and ask your ID department to install it. Good luck.
- Set output directories: The HPC scripts assume two locations for output storage:
    - One location for intermediate storage that can blow up quite a lot while scripts are running and also contains the raw sequencing reads – the main input for this pipeline. 
    - The `output/` directory of this repository for the permanent output, which will accumulate about 1.5 TB as you progress through the scripts.
- To set up the scripts accordingly, please adapt and run the following line:
```
intermediate="/absolute/path/to/folder" # WITHOUT SLASH AT THE END!
``` 
- Then navigate into the root folder of this repo and run:
```
repo_location=$(pwd)
sed -i "s|\$VSC_STAGING\/BPhage|${repo_location}|g" scripts/HPC/*.slrm
mkdir -p $repo_location/output/slurm_log

sed -i "s|\$VSC_SCRATCH\/BPhage|${intermediate}|g" scripts/HPC/*.slrm
mkdir -p $intermediate
```
- All set! I will refer to `$intermediate` in this README as the path to the intermediate storage but the variable does not have to remain set beyond this point.

### Download
- `download_raw_reads.slrm` (XXX): Download the raw data of this study from the SRA and rename the files. 
    - Requires: 
        - SRA accesstion list (XXX)
        - `data/rename.BPhage.nucleomics.sh` (XXX adapt for SRA numbers and move to `scripts/HPC`)
        - `data/rename.BPhage.nucleomics.scheme` (XXX adapt for SRA numbers)
    - Output: Nicely named fastq files with raw sequencing reads of all samples in `$intermediate/raw`.
- `download_additional_data.slrm`: Download the bee genome (also indexed here) and a bunch of additional data: 
    - Requires:
        - `data/bacteria.core.spec.tsv`
        - `data/other_datasets_SRA_accessions.tsv`
    - Output in `$intermediate/ref`: 
        - Nosema genome
        - Core bacteria genomes 
        - Apis mellifera genome (indexed with bowtie2)
    - Output in `$intermediate/additional_datasets`: 
        - Raw reads from other studies
        - Two INPHARED datasets 
        - Kirchberger et al's microvirus protines
    - More output: A clone of the MOP-UP github repository for the microvirus taxonomy downloaded into this repository (not tracked by git).

### Assembly
- `bphage_viper_with_dedup.slrm` (array of 471): ViPER pipeline
    - Requires: 
        - `data/BPhage.sample.list`
        - Raw read data in `$intermediate/raw`
    - Output: Trimmed reads, trimmed host-removed reads, assembly. For each sample there will be a separate folder in `$intermediate/bphage_viper_output`
- `reorganise_viper_output.slrm`: Re-organise ViPER output 
    - Requires: 
        - `data/BPhage.sample.list`
        - ViPER output in `$intermediate/bphage_viper_output`
    - Output: Oririnal files from each sample's viper output are moved into a common `CONTIGS`, `QC/FASTQC` (also with multiqc) `QC/QUAST` and `READ` (containing deduped trimmed and hostout reads) folder inside `output/bphage_viper_output`
- `scripts/HPC/cross_sample_clustering.slrm`: Cross-sample clustering 
    - Requires: Re-organised ViPER assemblies in `output/bphage_viper_output/CONTIGS/`
    - Output: Cross-sample clustered fasta and cluster member files: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`, `output/bphage_ALL_1kb_cross_95-85_clusters.tsv.gz`

### Phage identification
- `checkv.slrm`: CheckV 
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_checkv`
- `genomad.slrm`: geNomad
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_genomad`
- `filter_classification.slrm`: Filter for >= 50% complete phages: 
    - Requires: 
        - Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
        - CheckV output: `output/bphage_ALL_1kb_checkv/quality_summary.tsv.gz`
        - geNomad output: `output/bphage_ALL_1kb_genomad/bphage_ALL_1kb_cross_95-85_summary/bphage_ALL_1kb_cross_95-85_virus_summary.tsv.gz`
        - ICTV's VMR table: `data/VMR_19-250422_MSL37.xlsx` (downloaded from https://ictv.global/vmr on 14 March 2024)
        - Python script to do the actual filtering: `filter_classification.py` (called by the slrum script)
    - Output: 
        - Fasta file (`.fasta.gz`) and merged genomad and checkv table (`.csv`) with >= 50%-complete phages, >= 50%-complete unclassified viruses (containing contigs that are either "Unclassified" or classified only as "Viruses") and Picobirnaviridae (no completeness threshold):
            - `output/bphage_ALL_1kb_phages.*`, `output/bphage_ALL_1kb_unclassified_viruses.*`, `output/bphage_ALL_1kb_picobirna.*`

### Mapping
- `prepare_mapping.slrm`: Prepare mapping 
    - Requires: Fasta files from phage identification: `output/bphage_ALL_1kb_phages.fasta.gz`, `output/bphage_ALL_1kb_unclassified_viruses.fasta.gz`, `output/bphage_ALL_1kb_picobirna.fasta.gz`
    - Output: bwa-indexed mapping ref: `$intermediate/ref/bphage_mapping_ref.fasta`
- `bphage_mapping.slrm` (array of 471): Mapping  
    - Requires: 
        - `data/BPhage.sample.list`
        - bwa-indexed mapping ref: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: 
        - Mapping alignments and mapping stats of each sample in `output/mapped`
- `gather_mapping_stats.slrm`: Gather mapping stats 
    - Requires:
        - `data/BPhage.sample.list`
        - Mapping stats of individual samples in `output/mapped/`
        - Mapping ref to get the complete contig list: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: Separate tables for mapped reads, horizontal coverage and mean depth and a table with all contig lengths in `output/mapping_stats`

### Core contig refinement
- `contig_refinement.slrm` (array of 97): Prepare and run cobra 
    - Requires: 
        - List of core contigs: `data/core_contigs.txt`. These are defined later in the R part of the pipeline but were stored in `data` to avoid backtracking.
        - Individual assemblies pre-clustering: `output/bphage_viper_output/CONTIGS/*_all.contigs.fasta.gz`
        - Trimmed hostout reads from vipter ouput: `output/bphage_viper_output/READ/*.Hostout.R*.fastq.gz`
    - Output: Cobra-refined contigs at `output/core_contig_refinement/`
- `contig_refinement_stats.slrm`: Gather some stats 
    - Requires: Cobra outputs at `output/core_contig_refinement`
    - Output: 
        - Refinement stats: `output/core_contig_refinement/cobra_refinement_stats.tsv`
        - Fasta of extended contigs: `output/core_contig_refinement/extended_contigs.fasta`
- `contig_refinement_checkv_pharokka_phold.slrm`: Completeness and annotation
    - Requires:
        - Refinement stats: `output/core_contig_refinement/cobra_refinement_stats.tsv`
        - Fasta of extended contigs: `output/core_contig_refinement/extended_contigs.fasta`
        - Original CheckV output (for comparison): `output/bphage_ALL_1kb_checkv/quality_summary.tsv.gz`
        - Colabfold structures of proteins in `output/core_contig_refinement/colabfold_structures/extended_contigs_tophits/pdbs/` 
    - Output:
        - CheckV (incl. completeness comparison): `output/core_contig_refinement/extended_contigs_checkv`
        - Pharokka: `output/core_contig_refinement/extended_contigs_pharokka`
        - Phold: `output/core_contig_refinement/extended_contigs_phold`
        - Phold plots: `output/core_contig_refinement/extended_contigs_plots_phold`

### Annotation
- Pharokka: `scripts/HPC/annotation_pharokka_bphage.slrm`
    - Requires:
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
    - Output: `output/output/annotation/pharokka/bpgage_and_others`
- ColabFold: Done by collaborator George Bouras. Please refer to `scripts/HPC/protein_structures/README.md` for instructions.
    - Requires: Pharokka's prodigal-gv output: `output/annotation/pharokka/bpgage_and_others/bpgage_and_others_prodigal-gv.faa`
    - Output: Protein structures (AlphaFold2) of all identified proteins. Has to be placed into `output/annotation/phold_colabfold_structures`
- `annotation_phold_prepare.slrm`: Prepare phold compare
    - Requires: pdb files of ColabFold structures in `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
    - Output: 
        - Shortened pdb file names, otherwise phold will crash.
        - pdb files without corresponding entry in pharokka's output will be moved to `output/annotation/phold_colabfold_structures/basler_output_renamed/filtered_out_renamed_pdbs`
- `scripts/HPC/annotation_phold_compare.slrm`: Phold compare
    - Requires: 
        - Pharokka's gbk output: `output/annotation/pharokka_bphage_and_others/bphage_and_others.gbk`
        - Predicted structures at `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - sed script to re-lengthen contig names: `data/phold/relengthen.contig.names.sed`
        - Python script to filter gbk files (called from within the slrum script): `scripts/HPC/filter_gbk.py`
    - Output:
        - Predicted functions at: `output/annotation/phold_compare_bphage_and_others`
        - Additional files in this output folder with the original, long contig names: `output/annotation/phold_compare_bphage_and_others/*_long_names.*`
        - Ciros plots of all contigs, except a few very short ones that crashed `phold plot`: `output/annotation/plots_phold_compare_bphage_and_others` (using original, long contig names).
- `kegg_prepare.slrm`: Extract all "moron" gene sequences
    - Requires: Phold output: `output/core_contig_refinement/phold_compare_bphage_and_others/phold_aa_long_names.fasta`. Also from refined contigs: `output/core_contig_refinement/extended_contigs_phold/phold_aa_long_names.fasta`
    - Output: Gene sequences of "moron" genes: `output/kegg/moron.CDSs.fasta`. Use this to do a GhostKOALA search at the KEGG website (https://www.kegg.jp/ghostkoala/). The output has to be downloaded and turned into a useable format. As this is quite a manual work, the result of this process is at `data/kegg_mapping.tsv`

### Taxonomic classification
- `vcontact3_clustering_with_inphared.slrm`: vConTACT3
    - Requires: 
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
        - INPHARED dataset: `$intermediate/additional_datasets/inphared_3Apr2024_genomes_excluding_refseq.fa.gz`
    - Output: 
        - `output/vcontact3/bphage_vcontact3_b38_with_inphared/final_assignments.csv`
        - For visualisation: `output/vcontact3/bphage_vcontact3_b38_with_inphared/graph.bin_*.cyjs` (4 files). Run `scripts/R/cytoscape.R` on these files (see comments at the top of that script for instructions).
- `taxonomy_microviruses.slrm`: MOP-UP pipeline (for microviruses)
    - Requires: 
        - geNomad output: `output/bphage_ALL_1kb_genomad/bphage_ALL_1kb_cross_95-85_summary/bphage_ALL_1kb_cross_95-85_virus_summary.tsv.gz`
        - vConTACT3 output: `output/vcontact3/bphage_vcontact3_b38_with_inphared/final_assignments.csv`
        - Alternatively to the geNomad and vConTACT3 output: `/data/bphage.microviridae.contigs`
        - phold protein sequences: `output/annotation/phold_compare_bphage_and_others/phold_aa_long_names.fasta`
    - Output: `output/bphage_micros_mopup/bphage_micros_id30ForCytoscape.csv`. Run `scripts/R/microviruses_mopup_cytoscape.R` for visualisation (see comments at the top of this script for instructions).

### Lifestyle prediction
- `scripts/HPC/lifestyle_replidec.slrm`: Replidec 
    - Requires: 
        - Phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
        - Extended core contigs: `output/core_contig_refinement/extended_contigs.fasta`
- Output: Replidec's lifestyle prediction: `output/lifestyle/replidec/BC_predict.summary`

### Host prediction
- `iphop_bphage.slrm`: iPHoP
    - Requires: Phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
    - Output: Predicted hosts at: `output/host_prediction/`

### SNP analysis
- `SNP_analysis_mapping_prep.slrm`: Index bacterial genomes for mapping
    - Requires: Core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta`
    - Output: Indexed core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta.*`
- `SNP_analysis_mapping.slrm` (array of 150): Merge bam files per bee pool, retreive reads that didn't map against the phages and map those against the core bacteria genomes.
    - Requires: 
        - Bee pool list: `data/BPhage.bee.pool.list`
        - Indexed core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta.*`
        - Mapping output: `output/mapped_phages/`
    - Output: Reads mapping only to core bacteria genomes, not to phages or bee in (merged per bee pool) `$intermediate/merged_reads/`
- `SNP_analysis_SKA1.slrm` (array of 150): SKA fastq
    - Requires:
        - Bee pool list: `data/BPhage.bee.pool.list`
        - Trimmed reads (before host filtering): `output/bphage_viper_output/READ/*.trimmed.*.fastq.gz`
        - Hostout reads: `output/bphage_viper_output/READ/*.Hostout.*.fastq.gz`
    - Output: SKA kmer files for reads only mapping to bee, to phages or to core bacteria: `ska/skf_files/` (the the reads for the former two are deleted during cleanup)
- `SNP_analysis_distances_SKA1.slrm` (array of 3): SKA distance
    - Requires: SKA kmer files: `ska/skf_files/*.skf`
    - Output: Pairwise SNP distances in bee, phage and core bacteria read sets: `output/SNP_analysis/SKA_SNP_distances/*distances.tsv` 

### Nosema (Vairimorpha) mapping
- `Nosema_mapping_prep.slrm`: Index Nosema genomes for mapping
    - Requires: Nosema genome: `$intermediate/ref/Varimorpha_genomes.fasta`
    - Outoput: Indexed genome: `$intermediate/ref/Varimorpha_genomes.*`
- `Nosema_mapping_all.slrm` (array of 450): Mapping to Nosema genome
    - Requires: 
        - Sample list: `data/BPhage.sample.list`
        - Hostout reads: `output/bphage_viper_output/READ/*.Hostout.*.fastq.gz`
    - Output: 
        - Nosema-mapping alignments in `$intermediate/nosema_mapping_all`
        - Per-sample mapped read counts `$intermediate/nosema_mapping_all/*_read_counts.tsv`
- `Nosema_mapping_stats_all.slrm`: Gathering mapped read counts
    - Requires: Per-sample mapped read counts `$intermediate/nosema_mapping_all/*_read_counts.tsv`
    - Outout: Table of mapped reads: `output/nosema_mapped_counts_all.tsv`

### Additional datasets mapping
- `additional_datasets_mapping_with_unpaired.slrm` (array of 114): Mapping of reads from other studies to the phage genomes
    - Requires: 
        - List of SRA sccessions: `data/other_datasets_SRA_accessions.tsv`
        - Indexed phage genomes: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: Per-SRA coverage and mapped reads: `output/other_studies/*.coverage.gz`, `output/other_studies/*_read_stats.tsv`
- `additional_datasets_gather_mapping_stats_with_unpaired.slrm`: Gathering mapping stats
    - Requires:
        - Phage genomes: `$intermediate/ref/bphage_mapping_ref.fasta`
        - List of SRA sccessions: `data/other_datasets_SRA_accessions.tsv`
        - Per-SRA coverage and mapped reads: `output/other_studies/*.coverage.gz`, `output/other_studies/*_read_stats.tsv`
    - Outout: Stats for horizontal coverage, maped reads, mean depth and filtering stats: `output/other_studies/stats.other_studies.*`

### Novelty check: Clustering with INPHARED
- `inphared_clustering.slrm`: Clustering phage genomes with the INPHARED dataset on 70% identity over 85% of the genome length
    - Requires:
        - Phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
        - Inphared dataset: `$intermediate/additional_datasets/inphared_14Apr2025_genomes_excluding_refseq.fa.gz`
    - Output: Table with clusters: `output/inphared_clustering/bphage_and_inpha_70-85_clusters.tsv`

## R scripts
If you skipped the HPC part and jumped right here, I assume that you have cloned this repository and extracted the `mid_save.tar.gz` as descibed at the top of this README file.

If you worked through the HPC scripts, you will probably want another clone of this repository on a local computer and only carry over the output that is further needed. In that case, have a look at the contents of `mid_save.tar.gz` to see which files you will need:

```
tar -tf mid_save.tar.gz | grep -v "/$"
```

All the R scripts are meant to be run in RStudio in order of their numbering (for the scripts 8a, 8b etc. the order doesn't matter). While they will run start to finish without user interaction, I recommend you execute the commands one by one to trace potential issues. 

It would not be feasible to describe the in- and output of all scripts in detail here. Instead, I will provide general descriptions.

### R and package versions
- R 4.3.1
- RStudio 2023.12.1+402
- renv 1.1.4
- To reproduce R package versions run `renv::restore()`

### Filtering
`01.filtering.R`

Takes in the mapping stats in `output/mapping_stats_phages/` (mapped reads, horizontal coverage, mean depth) and sets mappped read counts of contigs to 0 if they have <70% horizontal coverage, <1 mean depth or are flagged as contaminants by `decontam`. It then writes updated read count tables and an updated geNomad classification table, filtering out samples and contigs that ended up with 0 reads after applying these filters.

---
### Alpha and beta diversity, relative abundances
`02.diversity_and_rel_abundance.R`

This is a very long script that takes about 2 hours to run. This is mainly due to the rarefaction (1.5 hours) and saving all the output plots (30 minutes). For a better overview, these parts are outsourced into separate scripts, which are called from within this script.

Alpha and beta diversities and relative abundances (called TPM in the scripts) are calculated and then many plots showing different aspects of the data are generated. The classification table is updated with vConTACT3 taxonomy and the CheckV output but the final version of this table is also stored as an R object in `data/classification.RDS`. The same is true for the metadata table (`data/metadata.RDS`). These R objects are loaded by downstream scripts whenever needed, so if you have different input, make sure to adapt those objects accordingly.

The following plots from this script appear as figures in the paper:
- **Figure 2**
    - `output/R/taxon_pies/pretty_pie.n.Class.*.pdf`
    - `output/R/taxon_pies/pretty_pie.n.Order.*.pdf`
    - `output/R/taxon_pies/pretty_pie.n.Family.*.pdf`
- **Figure 3a**
    - `output/R/prevalence/prevalence.Countries.pdf`
- **Figure 3b**
    - `output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/By_prevalence_Prevalence_Countries_relative_abundance.Country.pdf`
- **Supplementary Figure 1**
    - `output/R/beta/beta_all/Family_pcoa/beta.Family.all.all.pcoa.pdf`
    - `output/R/beta/beta_core_or_not/Family_pcoa/beta_core.no.Family.all.all.pcoa.pdf`
    - `output/R/beta/beta_core_or_not/Family_pcoa/beta_core.yes.Family.all.all.pcoa.pdf`
- **Supplementary Figure 4**
    - `output/R/absolute_counts_all_samples.pdf`
- **Supplementary Figure 5a**
    - `output/R/alpha/alpha_all/alpha_abs.Family.pdf`
    - `output/R/alpha/alpha_core_or_not/alpha_abs_core.no.Family.pdf`
    - `output/R/alpha/alpha_core_or_not/alpha_abs_core.yes.Family.pdf`
- **Supplementary Figure 5b**
    - `output/R/beta/beta_all/Family_pcoa/beta_abs.Family.all.all.pcoa.pdf`
    - `output/R/beta/beta_core_or_not/Family_pcoa/beta_abs_core.no.Family.all.all.pcoa.pdf`
    - `output/R/beta/beta_core_or_not/Family_pcoa/beta_abs_core.yes.Family.all.all.pcoa.pdf`
- **Supplementary Figure 8**
    - `output/R/rarefaction_thresholds/rarefaction_threshold.all.pdf`
    - `output/R/rarefaction_thresholds/rarefaction_threshold.noncore.pdf`
    - `output/R/rarefaction_thresholds/rarefaction_threshold.core.pdf`

---
`03.beta_dbRDA.R`

This script will pick up the the Bray-Curtis distance matrices and perform a distance-based redundancy analysis on them. This is done with the `RLdbRDA` package and cusomising their scripts to accept distance matrices instead of abundance tables.

The following plot from this script appear as figures in the paper:
- **Figure 3b**
    - `output/R/beta/beta_dbRDA/dbRDA.Family_patch.vertical.pdf`
---
`04_SNP_analysis.R` XXXX

---
### Lifestyle predictions
`05.lifestyle.R` XXXX

---
### Host predictions
`06.host.R`

This script takes the iPHoP output and makes plots out of it. Note that even though iPHoP was run with a confidence score threshold of 75, this script filters for confidence >90.

The following plots from this script appear as figures in the paper:
- **Figure 3d**
    - `output/R/host_pies/hosts.noncore.pdf`
    - `output/R/host_pies/hosts.core.pdf`
---
### Functional annotation
`07.gene_content.R`

This script combines the gene annotations from Phold, the KEGG assignments and the host predictions. Plots and values for PHROGs and KEGG assignments as well as metabolic gene prevalence are generated.

The following plots from this script appear as figures in the paper:
- **Figure 6a**
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.vertical.all.pdf`
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/legend.phrog.pdf`
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/kegg_bar.pdf`
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/legend.kegg.pdf`
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/goi_bar.pdf`
    - `output/R/genes_pathogens_and_landuse/phrog_and_kegg/legend.goi.pdf`
- **Figure 6b**
    - `output/R/genes_pathogens_and_landuse/hosts_of_genes_goi.pdf`
- **Supplementary Figure 2a**
    - `output/R/genes_pathogens_and_landuse/gene_prevalence.gene_facet.pdf`
- **Supplementary Figure 2b**
    - `output/R/genes_pathogens_and_landuse/gene_prevalence.overall.pdf`

---
### Mixed-effects models
`08a.mixed_models_gene_tpm_vs_landuse.R`, `08b.mixed_models_gene_presence_vs_landuse.R`, `08c.mixed_models_gene_tpm_vs_nosema_relabund.R`, `08d.mixed_models_pathogen_ct_vs_landuse.R`, `08e.mixed_models_pathogen_presence_vs_landuse.R`

These scripts are technically very similar. `a`, `c` and `d` perform linear mixed-effects models (LMMs) on relative gene abundances vs. land use parameters (`a`), on relative gene abundances vs. *Vairimorpha* relative abundance (`c`) and on pathogen Ct values vs. land use parameters (`d`). `b` and `e` perform generalized linear (logistic) mixed-effects models (GLMMs) on gene presence vs. land use parameters (`b`) and pathogen presence vs. land use parameters (`e`). 

The helper script `scripts/R/helpers/FAOstat_table.R` takes the country-wide pesticde usage (`data/FAOSTAT_pest_data_en_3-4-2025.csv`) and landuse (`data/FAOSTAT_area_data_en_3-5-2025.csv`) information from the FAO and combines it with the cropland area around the sampling sites (`data/land_cover_results.csv`) measured by COPERNICUS. Estimates of specific pesticide use at the sampling cites are then calculated. All numbers are from 2019, the year preceeding our sampling.

The landuse parameters all refer to a 2 km radius around the sampling sites. The cropland area was measured, pesticie uses are calculated estimates.
- Cropland area (ha)
    - Total pesticides
        - Insecticides
            - Insecticides – Organo-phosphates
            - Insecticides – Carbamates
            - Insecticides – Pyrethroids
            - Insecticides - nes
        - Herbicides
            - Herbicides – Phenoxy hormone products
            - Herbicides – Triazines
            - Herbicides – Amides
            - Herbicides – Carbamates
            - Herbicides – Dinitroanilines
            - Herbicides – Urea derivates
            - Herbicides - nes
        - Fungicides and Bactericides
            - Fung & Bact – Inorganics
            - Fung & Bact – Dithiocarbamates
            - Fung & Bact – Benzimidazoles
            - Fung & Bact – Triazoles, diazoles
            - Fung & Bact – Diazines, morpholines
            - Fung & Bact - nes
        - Plant Growth Regulators

In total there are 23 landuse parameters, 5 genes of interet (encoding PAPS reductasse, chitinase, glucosyltransferase, levanase and PnuC), 3 pathogens (BQCV, SBV and DWV-B) for the Ct value test in `d` and 3 pathogens (ABPV, V. ceranae and CBPV) for the pathogen presence/absence test in `e`. Benjamini-Hochberg correction of p-values was done in each script for all tests that successfully converged. Tests that failed to converge were excluded from further analyses. 

Script | Test | Number of tests | Successfully converge | Significant after BH correction
:---: | :---: | :---: | :---: | :---:
`a` | LMM (gene rel. abund vs. landuse) | 115 | 115 | 16
`b` | GLMM (gene presence vs. landuse)| 115 | 99 | 15
`c` | LMM (gene rel. abund vs. nosema rel. abund.)| 5 | 5 | 1
`d` | LMM (pathogen Ct vs. landuse) | 69 | 69 | 4
`e` | GLMM (paghogen presence vs. landuse) | 69 | 62 | 0

Several tables and plots are produced and placed into `output/R/genes_pathogens_and_landuse`, including plots of raw residuals vs. fitted values for model diagnostics. The figures in the paper were stiched together in the next script.

---
`09.mixed_models_main_figure.R` XXX

---
### Mapping to other datasets
`10_other_datasets.R` XXX

---
### Some additional calculations
`11_additional_small_tasks.R` XXX


The following plots from this script appear as figures in the paper:
- **Supplementary Figure 6**
    - `output/R/alpha/total_contigs_per_bee_pool.pdf`
- **Supplementary Figure 7**
    - `output/R/core_shared_between_guts.pdf`



