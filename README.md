# The honey bee triad: Phages are mutualistic partners in the gut microbiome of *Apis mellifera*

This repository contains all scripts and usage instructions to reproduce the analysis for the paper by Basler et al (XXX ref). 

This pipeline is split into two parts. The first part is meant for a high-performance computer (HPC) and can be skipped, if so wanted. There is also a reduced test dataset [available on Zenodo](https://zenodo.org/records/16937256) for the HPC part.

The second part is for the statistical analysis and visualisation using RStudio. If you want to skip the HPC part and only want to re-run the statistical analysis, please go straight to the ["R scripts"](#r-scripts) section of this README. The scripts from both parts pretend to be on the same computer but it is possible to clone this repo to an HPC and to a local computer, run the HPC scripts, copy the relevant output files from the HPC to the local computer and continue with the statistical analysis using the R project.

**Note**: The output of the tools and scripts will end up in the `output` folder inside the repo (which is why it's not tracked by git). The HPC scripts will create around 1.5 TB in total (or 1.5 GB if the test dataset is used), plus intermediate storage (see below). The R scripts create around 1 GB. Make sure to have enough free space or manage the output as it comes.

## Contents
- [HPC scripts](#hpc-scripts)
    - [General info](#general-info)
    - [Installations](#installations)
    - [Download](#download)
    - [Assembly](#assembly)
    - [Phage identification](#phage-identification)
    - [Mapping](#mapping)
    - [Core contig refinement](#core-contig-refinement)
    - [Annotation](#annotation)
    - [Taxonomic classification](#taxonomic-classification)
    - [Lifestyle prediction](#lifestyle-prediction)
    - [Host prediction](#host-prediction)
    - [SNP analysis](#snp-analysis)
    - [*Vairimorpha* (*Nosema*) mapping](#vairimorpha-nosema-mapping)
    - [Additional datasets mapping](#additional-datasets-mapping)
    - [Clustering with INPHARED](#clustering-with-inphared)
- [R scripts](#r-scripts)
    - [R and package versions](#r-and-package-versions)
    - [Package installations](#package-installations)
    - [Filtering](#filtering)
    - [Alpha and beta diversity, relative abundances](#alpha-and-beta-diversity-relative-abundances)
    - [Lifestyle predictions](#lifestyle-predictions)
    - [Host predictions](#host-predictions)
    - [Functional annotation](#functional-annotation)
    - [Associations with land use and pathogens (LMMs, GLMMs)](#associations-with-land-use-and-pathogens-lmms-glmms)
    - [Mapping to other datasets](#mapping-to-other-datasets)
    - [Accumulation curves](#accumulation-curves)
    - [Some additional calculations](#some-additional-calculations)
    - [Manuscript figures](#manuscript-figures)

## HPC scripts
### General info

- All scripts for this section are located in `scripts/HPC`.
- These scripts mostly represent jobs for a slurm scheduler. 
- You can run these scripts without using slrum. After following the installation instructions below, simply run the scripts with `bash -l <script.slrm>` instead of `sbatch <script.slrm>`. Please make sure you have the computing resources available that the scripts require (up to 36 cores and 144 GB RAM) or adjust the scripts accordingly. **Note:** If you want to avoid slrum, you will have to transform the array job scripts (meant to run in parallel) into successive loops. There is a script that can do that for you (see below), but this is only recommended for the test dataset.
- If you use slurm, you will have to adapt the instructions at the beginning of each script according to your setup. Particularly the `account` name and probably also the resource allocation will be different on your system.
- If you want to port these instructions to a different scheduler, your favourite generative AI tool will be very helpful.

### Installations
- To clone the repository, please run (v0.2.3 is frozen for review):
```
git clone --branch v0.2.3 --depth 1 https://github.com/nikolasbasler/BPhage
cd BPhage

```
- Mamba/conda: Please follow the official documentation (I recommend mamba: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). For this analysis, `mamba` v.1.4.2 was used.
- Conda environments: To install necessary conda environments, use the `.yml` files in `data/`. E.g. like this:
```
for env in data/env_*.yml; do 
    mamba env create -f $env
done

```
- ViPER pipeline: Please follow the instructions on Github **BUT** use the `viper_bphage` conda environment that you have just installed, instead of the environment from ViPER's Github page: https://github.com/Matthijnssenslab/ViPER. For this analysis, ViPER version 2.1 was used.
- For some tools, you will also need to install their databases. These are the databse versions used for this analysis:
 
 Tool | Database version | Source
 -- | -- | --
vConTACT3 | v223 | https://zenodo.org/records/10935513
Pharokka | v1.4.0 | https://zenodo.org/records/8267900
CheckV | v1.5 | https://portal.nersc.gov/CheckV
Phold | v0.2.0 | https://zenodo.org/records/12735568 
iPHOP | Aug_2023_pub_rw | use `iphop download --db_version iPHoP_db_Aug23_rw --db_dir /destination/path/`
geNomad | v1.7 | https://zenodo.org/records/10594875

- After downloading these databases, please adapt and run the following commands to provide the databases' locations:
```
vcontact_database="/absolute/path/to/database"
pharokka_database="/absolute/path/to/database"
checkv_database="/absolute/path/to/database"
phold_database="/absolute/path/to/database"
iphop_database="/absolute/path/to/database"
genomad_database="/absolute/path/to/database"

```

- **Note**: The tool MOP-UP (used for taxonomic clustering of the microviruses) relies on a library called `Boost`, which I didn't manage to install properly via conda but instead relied on a pre-installed module (installation happens in the `download_additional_data.slrm` script). If it gives you trouble, see version information in `data/env_mop-up_boost_module.txt` and ask your IT department to install it. If there is no way for you to install it, the rest of the pipeline still works. You will just not have Kirchberger et al's classification of microviruses. If you only run the test dataset, the microvirus classification can't be run anyway.
- Set output directories: The HPC scripts assume two locations for output storage:
    - One location for intermediate storage that can blow up quite a lot while scripts are running and also contains the raw sequencing reads – the main input for this pipeline. 
    - The `output/` directory of this repository for the permanent output, which will accumulate about 1.5 TB (or 1.5 GB if the test dataset is used) as you progress through the scripts.
- To set up the scripts accordingly, please adapt and run the following line (**without slash at the end!**):
```
intermediate="/absolute/path/to/folder"

``` 
- Then run the following commands to adapt the scripts to your folder paths. This will also incorporate the database paths given above.
```
repo_location=$(pwd)
sed -i "s|\$VSC_STAGING\/BPhage|${repo_location}|g" scripts/HPC/*.s*
mkdir -p $repo_location/output/slurm_log

sed -i "s|\$VSC_SCRATCH\/BPhage|${intermediate}|g" scripts/HPC/*.s*
mkdir -p $intermediate

sed -i "s|vcontact_db=.*|vcontact_db=\"${vcontact_database}\"|g" scripts/HPC/*.s*
sed -i "s|pharokka_db=.*|pharokka_db=\"${pharokka_database}\"|g" scripts/HPC/*.s*
sed -i "s|checkv_db=.*|checkv_db=\"${checkv_database}\"|g" scripts/HPC/*.s*
sed -i "s|phold_db=.*|phold_db=\"${phold_database}\"|g" scripts/HPC/*.s*
sed -i "s|iphop_db=.*|iphop_db=\"${iphop_database}\"|g" scripts/HPC/*.s*
sed -i "s|genomad_db=.*|genomad_db=\"${genomad_database}\"|g" scripts/HPC/*.s*

```
- If you want to turn the slrum array scripts into successive loops (only feasible with the test dataset) please also run the following. This will also reduce the number of datasets downloaded from other studies.
```
bash scripts/HPC/array_jobs_to_loop.sh

```
- If you did this, call these array scripts (or all scripts) with `bash -l <script.slrm>`, so the conda environments can get loaded from within those scripts.
- If running the test dataset, you also need to replace the sample and bee pool lists:
```
cat data/in_test_data > data/BPhage.sample.list 
cut -d "_" -f1-3 data/in_test_data > data/BPhage.bee.pool.list

```

- All set! I will refer to `$intermediate` in this README as the path to the intermediate storage but the variable does not have to remain set beyond this point.

### Download

- `download_raw_reads.slrm`: If you run the test dataset, skip this script and run `download_test_dataset.slrm` instead. Downloads the raw data of this study from the SRA and rename the files. There are 471 SRA datasets with 2 files each (forward and reverse reads). The total volume is >700 GB. I didn't parallelise this, because I don't think NCBI would allow hundreds of download requests from the same source. So expect this to take a while.
    - Requires: 
        - SRA accession list: `data/BPhage_SRAs.tsv`
    - Output: Nicely named fastq files with raw sequencing reads of all samples in `$intermediate/raw`.
- `download_test_dataset.slrm`: **Only** run if you want to use the test dataset!
    - Requries: Nothing
    - Output: Test dataset in `$intermediate/raw`. Note: If files from the full dataset are present, they will be overwritten!
- `download_additional_data.slrm`: Download the bee genome (also indexed here) and a bunch of additional data. If you have turned the array job scripts into successive loops (because you run the test dataset), only one SRA dataset will be downloaded here, otherwise mapping them later would take too much time.
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
- `bphage_viper_with_dedup.slrm` (array of 471): ViPER pipeline (read trimming, host removal, assembly). If you run the test dataset, skip this script and run `bphage_viper_with_dedup_TEST_DATA.slrm` instead.
    - Requires: 
        - Sample list: `data/BPhage.sample.list`
        - SRA accession list: `data/BPhage_SRAs.tsv`
        - Raw read data in `$intermediate/raw`
    - Output: Trimmed reads, trimmed host-removed reads, assembly. For each sample there will be a separate folder in `$intermediate/bphage_viper_output`
- `bphage_viper_with_dedup_TEST_DATA.slrm`: ViPER pipeline (see above).
- `reorganise_viper_output.slrm`: Re-organise ViPER output from one folder per sample to one folder per output.
    - Requires: 
        - `data/BPhage.sample.list`
        - ViPER output in `$intermediate/bphage_viper_output`
    - Output: Original files from each sample's viper output are moved into a common `CONTIGS`, `QC/FASTQC` (also with multiqc) `QC/QUAST` and `READ` (containing deduped trimmed and hostout reads) folder inside `output/bphage_viper_output`
- `cross_sample_clustering.slrm`: Collapsing redundancy by clustering contigs at 95% identity over 85% length of the smaller contig (see MIUVIG standard).
    - Requires: Re-organised ViPER assemblies in `output/bphage_viper_output/CONTIGS/`
    - Output: Cross-sample clustered fasta and cluster member files: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`, `output/bphage_ALL_1kb_cross_95-85_clusters.tsv.gz`

### Phage identification
- `checkv.slrm`: Genome completeness estimates with CheckV.
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_checkv`
- `genomad.slrm`: Virus identification with geNomad.
    - Requires: Cross-sample clustered contigs: `output/bphage_ALL_1kb_cross_95-85.fasta.gz`
    - Output: `output/bphage_ALL_1kb_genomad`
- `filter_classification.slrm`: Identify which virues are phages and filter for >= 50% complete phages.
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
- `prepare_mapping.slrm`: Index phage genomes for mapping.
    - Requires: Fasta files from phage identification: `output/bphage_ALL_1kb_phages.fasta.gz`, `output/bphage_ALL_1kb_unclassified_viruses.fasta.gz`, `output/bphage_ALL_1kb_picobirna.fasta.gz`
    - Output: bwa-indexed mapping ref: `$intermediate/ref/bphage_mapping_ref.fasta`
- `bphage_mapping.slrm` (array of 471): Mapping reads to the phage genomes. 
    - Requires: 
        - `data/BPhage.sample.list`
        - bwa-indexed mapping ref: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: 
        - Mapping alignments and mapping stats of each sample in `output/mapped`
- `gather_mapping_stats.slrm`: Gather mapping stats.
    - Requires:
        - `data/BPhage.sample.list`
        - Mapping stats of individual samples in `output/mapped/`
        - Mapping ref to get the complete contig list: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: Separate tables for mapped reads, horizontal coverage and mean depth and a table with all contig lengths in `output/mapping_stats`

### Core contig refinement
This cannot be done with the test dataset. If you are running the test dataset, please skip ahead to [Annotation](#annotation)

- `contig_refinement.slrm` (array of 97): Prepare and run cobra on core phage genomes to refine their assembly.
    - Requires: 
        - List of core contigs: `data/core_contigs.txt`. These are defined later in the R part of the pipeline but were stored in `data` to avoid backtracking.
        - Individual assemblies pre-clustering: `output/bphage_viper_output/CONTIGS/*_all.contigs.fasta.gz`
        - Trimmed hostout reads from vipter ouput: `output/bphage_viper_output/READ/*.Hostout.R*.fastq.gz`
    - Output: Cobra-refined contigs at `output/core_contig_refinement/`
- `contig_refinement_stats.slrm`: Gather some stats.
    - Requires: Cobra outputs at `output/core_contig_refinement`
    - Output: 
        - Refinement stats: `output/core_contig_refinement/cobra_refinement_stats.tsv`
        - Fasta of extended contigs: `output/core_contig_refinement/extended_contigs.fasta`
- `contig_refinement_checkv_pharokka_phold.slrm`: Completeness estimate and functional annotation of refined core phage genomes.
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
Only the first script (Pharokka) can be run with the test dataset. If you are running the test dataset, please skip ahead to [Taxonomic classification](#taxonomic-classification) after that.
- `annotation_pharokka_bphage.slrm`: First step of functional annotation with Pharokka.
    - Requires:
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
    - Output: `output/output/annotation/pharokka/bpgage_and_others`
- 3D structure prediction of all proteins with ColabFold: Done by George Bouras. Please refer to `scripts/HPC/protein_structures/README.md` for instructions.
    - Requires: Pharokka's prodigal-gv output: `output/annotation/pharokka/bpgage_and_others/bpgage_and_others_prodigal-gv.faa`
    - Output: Protein structures (AlphaFold2) of all identified proteins. Has to be placed into `output/annotation/phold_colabfold_structures`
- `annotation_phold_prepare.slrm`: Prepare phold compare.
    - Requires: pdb files of ColabFold structures in `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
    - Output: 
        - Shortened pdb file names, otherwise phold will crash.
        - pdb files without corresponding entry in pharokka's output will be moved to `output/annotation/phold_colabfold_structures/basler_output_renamed/filtered_out_renamed_pdbs`
- `annotation_phold_compare.slrm`: Look up predicted structures in a database with Phold compare.
    - Requires: 
        - Pharokka's gbk output: `output/annotation/pharokka_bphage_and_others/bphage_and_others.gbk`
        - Predicted structures at `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - sed script to re-lengthen contig names: `data/phold/relengthen.contig.names.sed`
        - Python script to filter gbk files (called from within the slrum script): `scripts/HPC/filter_gbk.py`
    - Output:
        - Predicted functions at: `output/annotation/phold_compare_bphage_and_others`
        - Additional files in this output folder with the original, long contig names: `output/annotation/phold_compare_bphage_and_others/*_long_names.*`
        - Ciros plots of all contigs, except a few very short ones that crashed `phold plot`: `output/annotation/plots_phold_compare_bphage_and_others` (using original, long contig names).
- `amg_phold_compare_with_temp_files.slrm`: Re-run phold on PAPS reductase genes and keep temp files to retreive more than just the top foldseek his.
    - Requires: 
        - Putative AMG CDSs: `data/AMG_CDSs.tsv`
        - Python script to filter down `.gbk` files: `grep_gbk.py` (called from within the script)
        - Pharokka output of non-refined contigs: `output/annotation/pharokka_bphage_and_others/bphage_and_others.gbk`
        - Pharokka output of refined contigs: `output/core_contig_refinement/extended_contigs_pharokka/extended_contigs.gbk`
        - Predicted structures at `output/annotation/phold_colabfold_structures/basler_output_renamed/renamed_pdbs/`
        - sed script to re-lengthen contig names: `data/phold/relengthen.contig.names.sed`
    - Output: Detailsed foldseek results: `output/annotation/phold_compare_bphage_and_others/paps_foldseek_results_long_names.tsv.gz`
- `kegg_prepare.slrm`: Extract all "moron" gene sequences.
    - Requires: Phold output: `output/core_contig_refinement/phold_compare_bphage_and_others/phold_aa_long_names.fasta`. Also from refined contigs: `output/core_contig_refinement/extended_contigs_phold/phold_aa_long_names.fasta`
    - Output: Gene sequences of "moron" genes: `output/kegg/moron.CDSs.fasta`. Use this to do a GhostKOALA search at the KEGG website (https://www.kegg.jp/ghostkoala/). The output has to be downloaded and turned into a useable format. As this is quite a manual work, the result of this process is at `data/kegg_mapping.tsv`

### Taxonomic classification
Only the first script (vConTACT3) can be run with the test dataset. If you are running the test dataset, please skip ahead to [Lifestyle predictions](#lifestyle-predictions) after that.

- `vcontact3_clustering_with_inphared.slrm`: Taxonomic classification of all genomes with vConTACT3. If you run the test dataset, I strongly recommend to skip this script and run `vcontact3_clustering_with_inphared_TEST_DATA.slrm` instead. It leaves out the INPHARED dataset and therefore runs much faster.
    - Requires: 
        - Phage, picobirna and unclassified bphage contigs: `output/bphage_ALL_1kb_phages.fasta.gz`, `bphage_ALL_1kb_picobirna.fasta.gz`, `bphage_ALL_1kb_unclassified_viruses.fasta.gz`
        - INPHARED dataset: `$intermediate/additional_datasets/inphared_3Apr2024_genomes_excluding_refseq.fa.gz`
    - Output: 
        - `output/vcontact3/bphage_vcontact3_b38_with_inphared/final_assignments.csv`
        - For visualisation: `output/vcontact3/bphage_vcontact3_b38_with_inphared/graph.bin_*.cyjs` (4 files). Run `scripts/R/cytoscape.R` on these files (see comments at the top of that script for instructions).
- `vcontact3_clustering_with_inphared_TEST_DATA.slrm`: vConTACT3 (see above)
- `taxonomy_microviruses.slrm`: MOP-UP pipeline for taxonomic classification of microviruses.
    - Requires: 
        - geNomad output: `output/bphage_ALL_1kb_genomad/bphage_ALL_1kb_cross_95-85_summary/bphage_ALL_1kb_cross_95-85_virus_summary.tsv.gz`
        - vConTACT3 output: `output/vcontact3/bphage_vcontact3_b38_with_inphared/final_assignments.csv`
        - Alternatively to the geNomad and vConTACT3 output: `/data/bphage.microviridae.contigs`
        - phold protein sequences: `output/annotation/phold_compare_bphage_and_others/phold_aa_long_names.fasta`
    - Output: `output/bphage_micros_mopup/bphage_micros_id30ForCytoscape.csv`. Run `scripts/R/microviruses_mopup_cytoscape.R` for visualisation (see comments at the top of this script for instructions).

### Lifestyle prediction
- `lifestyle_bacphlip`: Lifestyle prediction (virulent, temperate) with Bacphlip.
    - Requires: Phage, picobirna and unclassified contigs `$repo_location/output/bphage_ALL_1kb_phages_refined_contigs.fasta.gz`, `$repo_location/output/bphage_ALL_1kb_picobirna.fasta.gz`, `$repo_location/output/bphage_ALL_1kb_unclassified_viruses.fasta.gz`
    - Output: Bacphlip scores for virulent or temperate calls: `output/lifestyle/bphage_and_extended.fasta.bacphlip`

### Host prediction
- `iphop_bphage.slrm`: Host prediction with iPHoP.
    - Requires: Phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
    - Output: Predicted hosts at: `output/host_prediction/`

### SNP analysis
- `SNP_analysis_mapping_prep.slrm`: Index bacterial genomes for mapping.
    - Requires: Core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta`
    - Output: Indexed core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta.*`
- `SNP_analysis_mapping.slrm` (array of 150): Merge bam files per bee pool, retreive reads that didn't map against the phages (and therefore also not against the bee) and map those against the core bacteria genomes.
    - Requires: 
        - Bee pool list: `data/BPhage.bee.pool.list`
        - Indexed core bacteria genomes: `$intermediate/ref/core.bacteria.ref.genomes.fasta.*`
        - Mapping output: `output/mapped_phages/`
    - Output: Reads mapping only to core bacteria genomes, not to phages or bee in (merged per bee pool) `$intermediate/merged_reads/`
- `SNP_analysis_SKA1.slrm` (array of 150): Find slplit k-mers using SKA. If you run this script with the test dataset you will get errors because of missing files. This is expected and can be ignored.
    - Requires:
        - Bee pool list: `data/BPhage.bee.pool.list`
        - Trimmed reads (before host filtering): `output/bphage_viper_output/READ/*.trimmed.*.fastq.gz`
        - Hostout reads: `output/bphage_viper_output/READ/*.Hostout.*.fastq.gz`
    - Output: SKA kmer files for reads only mapping to bee, to phages or to core bacteria: `ska/skf_files/` (the the reads for the former two are deleted during cleanup)
- `SNP_analysis_distances_SKA1.slrm` (array of 3): Calculate pair-wise SNP distances between samples using SKA distance. This is done independently for phage reads, bacteria reads and bee reads.
    - Requires: SKA kmer files: `ska/skf_files/*.skf`
    - Output: Pairwise SNP distances in bee, phage and core bacteria read sets: `output/SNP_analysis/SKA_SNP_distances/*distances.tsv` 

### *Vairimorpha* (*Nosema*) mapping
- `Nosema_mapping_prep.slrm`: Concatenate and index Nosema genomes for mapping.
    - Requires: Nosema genome: `$intermediate/ref/Varimorpha_genomes.fasta`
    - Outoput: Indexed genome: `$intermediate/ref/Varimorpha_genomes.*`
- `Nosema_mapping_all.slrm` (array of 450): Mapping reads to Nosema genome.
    - Requires: 
        - Sample list: `data/BPhage.sample.list`
        - Hostout reads: `output/bphage_viper_output/READ/*.Hostout.*.fastq.gz`
    - Output: 
        - Nosema-mapping alignments in `$intermediate/nosema_mapping_all`
        - Per-sample mapped read counts `$intermediate/nosema_mapping_all/*_read_counts.tsv`
- `Nosema_mapping_stats_all.slrm`: Gathering mapped read counts of nosema mapping.
    - Requires: Per-sample mapped read counts `$intermediate/nosema_mapping_all/*_read_counts.tsv`
    - Outout: Table of mapped reads: `output/nosema_mapped_counts_all.tsv`

### Additional datasets mapping
- `additional_datasets_mapping_with_unpaired.slrm` (array of 114): Mapping of reads from other studies to the phage genomes. If you turned this array job script into a successive loop, only one of the datasets are being mapped, otherwise it would take too long, even with the test dataset.
    - Requires: 
        - List of SRA sccessions: `data/other_datasets_SRA_accessions.tsv`
        - Indexed phage genomes: `$intermediate/ref/bphage_mapping_ref.fasta`
    - Output: Per-SRA coverage and mapped reads: `output/other_studies/*.coverage.gz`, `output/other_studies/*_read_stats.tsv`
- `additional_datasets_gather_mapping_stats_with_unpaired.slrm`: Gathering mapping stats.
    - Requires:
        - Phage genomes: `$intermediate/ref/bphage_mapping_ref.fasta`
        - List of SRA sccessions: `data/other_datasets_SRA_accessions.tsv`
        - Per-SRA coverage and mapped reads: `output/other_studies/*.coverage.gz`, `output/other_studies/*_read_stats.tsv`
    - Outout: Stats for horizontal coverage, maped reads, mean depth and filtering stats: `output/other_studies/stats.other_studies.*`

### Clustering with INPHARED
- `inphared_clustering.slrm`: Clustering phage genomes with the INPHARED dataset on 70% identity over 85% of the genome length.
    - Requires:
        - Phage, picobirna and unclassified contigs at `output/bphage_ALL_1kb_*.fasta.gz`
        - Inphared dataset: `$intermediate/additional_datasets/inphared_14Apr2025_genomes_excluding_refseq.fa.gz`
    - Output: Table with clusters: `output/inphared_clustering/bphage_and_inpha_70-85_clusters.tsv`

### Metabolic gene confirmation
- `amg_confirmation_get_annotations.slrm`
    - Requires:
        - List of genomes with potential AMGs: `data/AMG_genomes.tsv`
        - Pharokka's output of the bulk contigs: `output/annotation/pharokka_bphage_and_others/prodigal-gv.faa`
        - Pharokka's output of the refined contigs: `output/core_contig_refinement/extended_contigs_pharokka/prodigal-gv.faa`
        - Pyhon script to filter `.gbk` files (will be called from within the script): `scripts/HPC/grep_gbk.py` 
        - List of CDSs of potential AMGs: `data/AMG_CDSs.tsv`
    - Output: Genome sequences of contigs with potential AMGs and their protein sequences in `output/amg_confirmation`

---
---

## R scripts
If you skipped the HPC part and jumped right here, you will want to clone this repository to a computer that runs RStudio and then extracte the `mid_save.tar.gz`. Extracting this file with `-k` will not overwrite existing files, so it's safe to use if you generated some HPC output. v0.2.3 is frozen for review.

```
git clone --branch v0.2.3 --depth 1 https://github.com/nikolasbasler/BPhage
cd BPhage
tar -kxvzf mid_save.tar.gz

```

If you worked through the HPC scripts, you will probably want another clone of this repository on a local computer and only carry over the output that is further needed. In that case, have a look at the contents of `mid_save.tar.gz` to see which files you will need:

```
tar -tf mid_save.tar.gz | grep -v "/$"

```

All the R scripts are meant to be run in RStudio in order of their numbering. Each script can run start to finish without user interaction in a few seconds, except `02.diversity_and_rel_abundance.R`, which takes about 1h and `03.beta_dbRDA.R`, which takes about 15 minutes. I recommend to restart the RStudio session before every script.

It would not be feasible to describe the in- and output of all scripts in detail here. Instead, I will provide general descriptions. [At the end of this README](#manuscript-figures), there is a table with all figures that appear in the paper, the script that creates them and their file names.

### R and package versions
- R 4.3.1
- RStudio 2023.12.1+402
- renv 1.1.4

### Package installations
- For reproducibility, it would be best to install R 4.3.1. On Windows, RStudio supports switching between different R versions. On Mac or Linux, you may have to rely on `rig` (https://github.com/r-lib/rig) to do that.
- Once you open the R project (`BPhage.Rproj`), `renv` should install itself. If not, please install it yourself.
- To reproduce R package versions run `renv::restore()` in the RStudio terminal.
    - Don't worry about the `ERROR [error code 22]` messages during download. `renv` will keep trying and is usually able to download a package after a few attempts.
    - The message `GitHub authentication credentials are not available` can also be ignored.
- Restart the RStudio session after successful installation
- If the pacakge installation fails, you may have to install a missing compiler. Unfortunately, I can't provide instructions for all possible issues here, so you will have to resolve these yourself. Your trusted AI chatbot might be of service here. Alternatively, you can of course install the packages manually, but there are many and if the versions don't line up with the ones used here, the scripts might crash.

### Filtering
`01.filtering.R`

Takes in the mapping stats in `output/mapping_stats_phages/` (mapped reads, horizontal coverage, mean depth) and sets mappped read counts of contigs to 0 if they have <70% horizontal coverage, <1 mean depth or are flagged as contaminants by `decontam`. It then writes updated read count tables and an updated geNomad classification table, filtering out samples and contigs that ended up with 0 reads after applying these filters.

---
### Alpha and beta diversity, relative abundances
`02.diversity_and_rel_abundance.R`

This is a very long script that takes about 1 hour to run. This is mainly due to the rarefaction (45 minutes) and saving all the output plots (10 minutes). For a better overview, these parts are outsourced into separate scripts, which are called from within this script. For the rarefaction (1000 iterations), I make use of the `future` package to run on 4 CPU cores in parallel. If you want to change this, adapt the `plan(multisession, workers=4)` calls in the helper scripts `scripts/R/helpers/alpha_diversity.R` and `scripts/R/helpers/beta_diversity.R`.

Alpha and beta diversities and relative abundances (called TPM in the scripts) are calculated and then many plots showing different aspects of the data are generated. This mainly served the purpose of data exploration. The `classification` table is updated with `vConTACT3` taxonomy and the `CheckV` output but the final version of the `classification` table is also stored as an R object in `data/classification.RDS`. The same is true for the metadata table (`data/metadata.RDS`). These R objects are loaded by all scripts whenever needed, so if you have different input, make sure to adapt those objects accordingly.

---

`02b.abs_rel_diversity_correlations.R`

Making correlations and repective plots between Shannon and Bray-Curtis diversities from read data vs. from abolute count data. 

---


`03.beta_dbRDA.R`

Picks up the Bray-Curtis dissimilarity matrices and performs distance-based redundancy analyses on them. Runs for about 15 minutes.

---
`04_SNP_analysis.R`

PCoA and distance-based redundancy analysis is performed on SKA's pairwise SNP distances, followed by Mantel tests.

---
### Lifestyle predictions
`05.lifestyle.R`

Replidec's lifestyle predictions are summed up in this script. The `classification` table is updated accordingly.

---
### Host predictions
`06.host.R`

This script takes the iPHoP output and makes plots out of it. Note that even though iPHoP was run with a confidence score threshold of 75, this script filters for confidence >90. The `classification` table is updated accordingly.

---
### Functional annotation
`07.gene_content.R`

This script combines the gene annotations from Phold, the KEGG assignments and the host predictions. Plots and values for PHROGs and KEGG assignments as well as metabolic gene prevalence are generated.

---
### Associations with land use and pathogens (LMMs, GLMMs)
`08a.mixed_models_gene_tpm_vs_landuse.R`, `08b.mixed_models_gene_presence_vs_landuse.R`, `08c.mixed_models_pathogen_ct_vs_landuse.R`, `08d.mixed_models_pathogen_presence_vs_landuse.R`

These scripts are technically very similar. `a` and `c` perform linear mixed-effects models (LMMs) on relative gene abundances vs. land use parameters (`a`), and on pathogen Ct values vs. land use parameters (`c`). `b` and `d` perform generalized linear (logistic) mixed-effects models (GLMMs) on gene presence vs. land use parameters (`b`) and pathogen presence vs. land use parameters (`d`). 

The helper script `scripts/R/helpers/FAOstat_table.R` takes the country-wide pesticde usage (`data/FAOSTAT_pest_data_en_3-4-2025.csv`) and landuse (`data/FAOSTAT_area_data_en_3-5-2025.csv`) information from the FAO and combines it with the cropland area around the sampling sites (`data/land_cover_results.csv`) measured by COPERNICUS. Estimates of specific pesticide use at the sampling cites are then calculated. All numbers are from 2019, the year preceeding our sampling.

The landuse parameters all refer to a 2 km radius around the sampling sites. The cropland area was measured, pesticie uses are calculated estimates.
1. Cropland area
1. Total pesticides (subsumes all pesticide use)
1. Insecticides (use of all insecticide sub-groups)
1. Insecticides – Organo-phosphates
1. Insecticides – Carbamates
1. Insecticides – Pyrethroids
1. Insecticides - nes
1. Herbicides (use of all herbicide sub-groups)
1. Herbicides – Phenoxy hormone products
1. Herbicides – Triazines
1. Herbicides – Amides
1. Herbicides – Carbamates
1. Herbicides – Dinitroanilines
1. Herbicides – Urea derivates
1. Herbicides - nes
1. Fungicides and Bactericides (use of all fungicide & bactericide sub-groups)
1. Fung & Bact – Inorganics
1. Fung & Bact – Dithiocarbamates
1. Fung & Bact – Benzimidazoles
1. Fung & Bact – Triazoles, diazoles
1. Fung & Bact – Diazines, morpholines
1. Fung & Bact - nes
1. Plant Growth Regulators

In total there are 23 landuse parameters, 1 genes of interet (encoding PAPS reductasse), 3 pathogens (BQCV, SBV and DWV-B) for the Ct value test in `c` and 3 pathogens (ABPV, V. ceranae and CBPV) for the pathogen presence/absence test in `d`. Benjamini-Hochberg correction of p-values was done in each script for all tests that successfully converged. Tests that failed to converge were excluded from further analyses. 

Script | Test | Number of tests | Successfully converged | Significant after BH correction
:---: | :---: | :---: | :---: | :---:
`a` | LMM (gene rel. abund. vs. landuse) | 23 | 23 | 13
`b` | GLMM (gene presence vs. landuse)| 23 | 21 | 3
`c` | LMM (pathogen Ct vs. landuse) | 69 | 69 | 4
`d` | GLMM (paghogen presence vs. landuse) | 69 | 62 | 0

Several tables and plots are produced and placed into `output/R/genes_pathogens_and_landuse` (and subfolders), including plots of raw residuals vs. fitted values for model diagnostics. The figures in the paper were stiched together in the next script.

---
`09.mixed_models_figures.R`
This script takes the figures from the mixed model scripts and puts them together in different variations.

---
### Mapping to other datasets
`10_other_datasets.R`

From the mapping of the other studies' SRA datasets, the same filters are applied as for the mapping of this study. The mapping stats are then taken to make several tables and an upset plot showing presence of core phages in those studies. 

---
### Accumulation curves
`11_accumulation_curves.R`

Several accumulation curves of phage genomes over all samples, bee pools or specific gut parts are generated.

---
### Some additional calculations
`12_additional_small_tasks.R`

Some smaller tasks that don't deserve their own scripts are done here. For example, the `classification` table is updated with the INPHARED clustering information.

### Manuscript figures

Figure | Created in script | Files (inside `output/R/`) 
--- | --- |--- 
Figure 1 | NA | NA
Figure 2 | `02.diversity_and_rel_abundance.R` | `taxon_pies/pretty_pie.n.Class.pie.pdf` <br> `taxon_pies/pretty_pie.n.Class.legend.pdf` <br> `taxon_pies/pretty_pie.n.Order.pie.pdf` <br> `taxon_pies/pretty_pie.n.Order.legend.pdf` <br> `taxon_pies/pretty_pie.n.Family.pie.pdf` <br> `taxon_pies/pretty_pie.n.Family.legend.pdf`
Figure 3a | `02.diversity_and_rel_abundance.R` | `prevalence/prevalence.Countries.pdf`
Figure 3b | `02.diversity_and_rel_abundance.R` | `relative_abundance/relative_abundance_by_metavar_core_or_not/By_prevalence_Prevalence_Countries_relative_abundance.Country.pdf`
Figure 3c | `05.lifestyle.R` | `lifestyle/replidec.Caudoviricetes.all.horizontal.pdf`
Figure 3d | `06.host.R` | `host_pies/hosts.noncore.pdf` <br> `host_pies/hosts.core.pdf`
Figure 4a | NA | NA
Figure 4b | `10_other_datasets.R` | `other_studies/core_read_presence_overlap.upset.patch.pdf`
Figure 5a | `02.diversity_and_rel_abundance.R` | `alpha/pretty_alpha_selection.pdf`
Figure 5b | `03.beta_dbRDA.R` | `beta/beta_dbRDA/dbRDA.Family_patch.vertical.pdf`
Figure 5c | `04_SNP_analysis.R` | `SNP_analysis/SNP_RDA_horizontal.pdf`
Figure 6a | `07.gene_content.R` | `genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.vertical.all.pdf` <br> `genes_pathogens_and_landuse/phrog_and_kegg/legend.phrog.pdf` <br> `genes_pathogens_and_landuse/phrog_and_kegg/kegg_bar.pdf` <br> `genes_pathogens_and_landuse/phrog_and_kegg/legend.kegg.pdf` <br> `genes_pathogens_and_landuse/phrog_and_kegg/goi_bar.pdf` <br> `genes_pathogens_and_landuse/phrog_and_kegg/legend.goi.pdf`
Figure 6b | `07.gene_content.R` | `genes_pathogens_and_landuse/hosts_of_genes_goi.pdf`
Figure 6c | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/paps_tpm.pdf`
Figure 6d | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/goi_presence.pdf`
Figure 6e | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/nosema_relabund_and_ct.pdf`
Supplementary Figure 1 | `02.diversity_and_rel_abundance.R` | `beta/beta_all/Family_pcoa/beta.Family.all.all.pcoa.pdf` <br> `beta/beta_core_or_not/Family_pcoa/beta_core.no.Family.all.all.pcoa.pdf` <br> `beta/beta_core_or_not/Family_pcoa/beta_core.yes.Family.all.all.pcoa.pdf`
Supplementary Figure 2a | `07.gene_content.R` | `genes_pathogens_and_landuse/gene_prevalence.gene_facet.pdf`
Supplementary Figure 2b | `07.gene_content.R` | `genes_pathogens_and_landuse/gene_prevalence.overall.pdf`
Supplementary Figure 3a | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/sup_fig_presence.pdf` <br> `genes_pathogens_and_landuse/selected_graphs/sup_fig_presence_legend.pdf`
Supplementary Figure 3b | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/sup_fig_tpm.pdf` <br> `genes_pathogens_and_landuse/selected_graphs/sup_fig_tpm_legend.pdf`
Supplementary Figure 3c | `09.mixed_models_figures.R` | `genes_pathogens_and_landuse/selected_graphs/sup_fig_ct.pdf` <br> `genes_pathogens_and_landuse/selected_graphs/sup_fig_ct_legend.pdf`
Supplementary Figure 4 | `02.diversity_and_rel_abundance.R` | `absolute_counts_all_samples.pdf`
Supplementary Figure 5a | `02.diversity_and_rel_abundance.R` | `alpha/pretty_alpha_selection_abs.pdf` 
Supplementary Figure 5b | `02b.abs_rel_diversity_correlations.R` | `alpha/rel_abs_shannon_correlation/alpha_rel_abs_shannon_correlation.wrap.pdf`
Supplementary Figure 5c | `02.diversity_and_rel_abundance.R` | `beta/beta_all/Family_pcoa/beta_abs.Family.all.all.pcoa.pdf` <br> `beta/beta_core_or_not/Family_pcoa/beta_abs_core.no.Family.all.all.pcoa.pdf` <br> `beta/beta_core_or_not/Family_pcoa/beta_abs_core.yes.Family.all.all.pcoa.pdf`
Supplementary Figure 5d | `02b.abs_rel_diversity_correlations.R` | `beta/rel_abs_bray_correlation/beta_rel_abs_bray_correlation.wrap.pdf`
Supplementary Figure 6 | `11_accumulation_curves.R` | `accumulation_curves/accumulation_curve.all.all.pdf`
Supplementary Figure 7 | `12_additional_small_tasks.R` | `core_shared_between_guts.pdf`
Supplementary Figure 8 | `02.diversity_and_rel_abundance.R` | `rarefaction_thresholds/rarefaction_threshold.all.pdf` <br> `rarefaction_thresholds/rarefaction_threshold.noncore.pdf` <br> `rarefaction_thresholds/rarefaction_threshold.core.pdf`


