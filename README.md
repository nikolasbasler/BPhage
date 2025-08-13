# The honey bee triad: Phages are mutualistic partners in the gut microbiome of *Apis mellifera*

This repository contains all scripts to reproduce the analysis for the paper by Basler et al (XXX ref) and usage instructions. The entire 2,343 assembled phage genomes are also stored here as well as a subset with the 97 core phages. The 1,066 >90% complete genomes are on Genbank (See Supplementary file, sheet i - Genbank accessions).

This pipeline is split into two parts. The first part is meant for a high-performance computer (HPC) and can be skipped, if so wanted. The second part is for the statistical analysis and visualisation using RStudio. The idea of the workflow is to set up on the HPC, run the respective scripts, clone this repo to a local computer, copy the relevant output files from the HPC to the local computer and continue with the statistical analysis using the R project. 

To clone the repository, please run
```
git clone https://github.com/nikolasbasler/BPhage
```
**Note**: The output of the tools and scripts will end up in the `output` folder inside the repo (which is why it's not tracked by git). The HPC scripts will create around 1.5 TB in total (plus intermediate storage, see below), the R scripts around 1 GB. Make sure to have enough free space or manage the output as it comes.

**If you want to skip the HPC part and only want to re-run the statistical analysis:** Clone the repo as mentioned above, extract the mid-save archive, and then skip ahead to the "R scripts" section of this readme file.
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

### Installations (manual, no scripts)
For tool versions see the conda environment `.yml` files: `data/env_*.yml`
- Mamba/conda: Please follow the official documentation (e.g.: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). For this analysis, `mamba` v.1.4.2 and `conda` v.23.1.0 were used.
- ViPER pipeline: Please follow the instructions on Github **BUT** use `data/env_viper_bphage.yml` from this repository instead of ViPER's `viper.yml`: https://github.com/Matthijnssenslab/ViPER. For this analysis, ViPER version 2.1 was used.
- Additional conda environments: To install additional conda environments with all necessary tools, use the `env_*.yml` files in `data/`: `env_cobra.yml`, `env_iphop.yml`, `env_mopup.yml`, `env_ncbi.yml`, `env_pharokka.yml`, `env_phold.yml`, `env_replidec.yml`, `env_ska.yml`, `env_vcontact3.yml`. E.g. like this:
```
mamba env create -f env_cobra.yml
```

- **Note**: The tool MOP-UP (used for taxonomic clustering of the microviruses) relies on a library called `Boost`, which I didn't manage to install properly via conda but instead relied on a pre-installed module. If it gives you trouble, see version information in `env_mop-up_boost_module.txt`.
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
- All set! I will refer to `$intermediate` in this readme as the path to the intermediate storage but the variable does not have to remain set beyond this point.

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

### Annotation
- Pharokka: `scripts/HPC/annotation_pharokka_bphage.slrm`
    - Requires:
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
    - Output: `output/output/annotation/pharokka/bpgage_and_others`

#### CONTINUE HERE. WAIT FOR PHAROKKA TO FINISH ON THE SET WITHOUT THE OTHER DATASET'S CONTIGS

- Collabfold: Done by collaborator George Bouras.
    - Requires: Pharokka's prodigal-gv output: `output/annotation/pharokka/bpgage_and_others/prodigal-gv.faa` (in communication with George the file was named `bpgage_and_others_prodigal-gv.faa`)
    - Output: Protein structures (Alpha fold 2) of all identified proteins from bphage and other studies. Provided by George, placed into `output/annotation/phold_colabfold_structures`
- `annotation_phold_prepare.slrm`: Prepare phold compare
    - Requires: 
        - pdb files of George's collabfold structures in `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - List of contigs that were sent to George but are not in the bphage_and_others dataset anymore because they were either filtered out by geNomad/CheckV (from additional studies) or collapsed by clustering: `$repo_location/data/phold/sent_to_george_but_removed_from_original_dataset`
    - Output: 
        - Renamed pdb files with shortened names, otherwise phold will crash.
        - pdb files without corresponding entry in pharokka's output will be moved to `output/annotation/phold_colabfold_structures/basler_output_renamed/filtered_out_renamed_pdbs`
- Phold: `scripts/HPC/annotation_phold_compare.slrm`
    - Requires: 
        - Pharokka's gbk output: `output/annotation/pharokka_bphage_and_others/bphage_and_others.gbk`
        - Predicted structures at `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - sed script to re-lengthen contig names: `data/phold/relengthen.contig.names.sed`
    - Output:
        - Predicted functions at: `output/annotation/phold_compare_bphage_and_others`
        - Additional files in this output folder with the original, long contig names: `output/annotation/phold_compare_bphage_and_others/*_long_names.*`
        - Ciros plots of all contigs, except a few very short ones that crashed `phold plot`: `output/annotation/plots_phold_compare_bphage_and_others` (using original, long contig names).

### Taxonomic clustering
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

### Lifestyle prediction
- `scripts/HPC/lifestyle_replidec.slrm`: Replidec 
    - Requirements: BPhage phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
    - Extended core contigs: `output/core_contig_refinement/extended_contigs.fasta`
- Output: Replidec's lifestyle prediction: `output/lifestyle/replidec/BC_predict.summary`

### Host prediction


### SNP analysis
- mapping prep (bact)
- mapping vs bact
- SKA `SNP_analysis_SKA1`
- SKA distances
<!-- - mapping prep
- mapping (array of 150)
- (subspecies trimming)
- SKA ()
- (SKA bacteria)
- (SKA subspecies)
- SKA distances -->

### Nosema mapping

### Additional datasets mapping

### Novelty check: Clustering with INPHARED

## R scripts
### Package versions
- R 4.3.1
- RStudio 2023.12.1+402
- renv 1.1.4
- To reproduce R package versions run `renv::restore()`

- Filtering script `xxx`
    - Requires: 
        - `xxx`
    - Output:
        - `xxx`
        - Filtered classifiaction table: `output/R/phage.filt.gnmd.classification.csv`
- Visualise predictions and define host groups


