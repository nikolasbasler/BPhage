# Do manually by copy-and-pasting line by line
# This is not part of the pipeline but included here in case it needs to be remembered. 
# These issues may be fixed with future versions of Pharokka, Phold and suvtk (https://github.com/LanderDC/suvtk)

repo_location="$VSC_STAGING/BPhage"

cd $VSC_STAGING/submission_genbank_with_core/

conda activate suvtk_github

# suvtk download-database

# Shorten the names
# Also remove lines with the following patterns:
#   host specificity=
#   inhibitor of host transcription=
#   inhibits RNA polymerase=
#   negatively on RNaseE=
#   positively regulates RNaseIII=
#   target for F exclusion=
sed -E '/NODE_.*_length_.*_cov_/s/(NODE_|length_|cov_)//g' long_names/taxonomy.tsv > s_taxonomy.tsv
sed -E '/NODE_.*_length_.*_cov_/s/(NODE_|length_|cov_)//g' long_names/checkv.tsv > s_checkv.tsv
sed -E '/NODE_.*_length_.*_cov_/s/(NODE_|length_|cov_)//g' long_names/sourcefile.tsv > s_sourcefile.tsv 
sed -E '/NODE_.*_length_.*_cov_/s/(NODE_|length_|cov_)//g' long_names/genomes.list > s_genomes.list
sed -E '/NODE_.*_length_.*_cov_/s/(NODE_|length_|cov_)//g' $repo_location/output/annotation/phold_compare_bphage_and_others/phold_long_names.gbk | \
    grep -v -e "host specificity=" -e "inhibitor of host transcription=" -e "inhibits RNA polymerase=" -e "negatively on RNaseE=" -e "positively regulates RNaseIII=" -e "target for F exclusion=" \
    > s_all_genomes.gbk

# Remove locus tag (the tag is not really needed and still gives warnings for too long names)
python3 filter_gbk_and_remove_locus_tag.py s_all_genomes.gbk s_genomes.list s_genomes.gbk

suvtk virus-info --taxonomy s_taxonomy.tsv --database /staging/leuven/stg_00029/DB/suvtk_db --output virus_info

mkdir -p genbank_to_table
suvtk gbk2tbl -i s_genomes.gbk -p genbank_to_table/genbank_to_table

suvtk comments -t virus_info/miuvig_taxonomy.tsv -f for_comments_call/features.tsv \
    -m for_comments_call/global_miuvig.tsv -a for_comments_call/assembly.tsv \
    -o structured_comment -c s_checkv.tsv

# Wrong qualifiers or values:
#   pseudogene
#   transl_table for tRNA and rmRNA entries
#   source for Pyrodigal should be inferance with ab initio prediction. Also version number should be separated with : instead of _
#   source for tRNAscan should be inference with profile. Also version number should be separated with : instead of _
#   source for Aragorn should be inference with profile. Also version number should be separated with : instead of _
#   source for minced should be inference with nucleotide motif. Here the version number is separated with :
sed 's/pseudogene/tRNA\n\t\t\tpseudogene\tunknown/g' genbank_to_table/genbank_to_table.tbl |
        awk -F'\t' '
        # 1) If flag is set and this line contains "transl_table", delete it and clear flag
        flag && /transl_table/ { flag=0; next }

        # 2) If this is a tRNA/tmRNA row, set flag (and print it)
        $3=="tRNA" || $3=="tmRNA" || $3=="repeat_region" { flag=1; print; next }

        # 3) If the 3rd column of this row is something else, clear flag
        $3!="" && $3!="tRNA" && $3!="tmRNA" && $3!="repeat_region" { flag=0 }

        # 4) Otherwise just print
        { print }
        ' | \
        sed 's/source\tPyrodigal-gv_0.3.1/inference\tab initio prediction:Pyrodigal-gv:0.3.1/g' | \
        sed 's/source\ttRNAscan-SE_2.0.12/inference\tprofile:tRNAscan-SE:2.0.12/g' | \
        sed 's/source\tAragorn_1.2.41/inference\tprofile:Aragorn:1.2.41/g' | \
        sed 's/source\tminced:0.4.2/inference\tnucleotide motif:minced:0.4.2/g' \
        > genbank_to_table/genbank_to_table_fixed.tbl 

# Run this to catch errors and parse through them
suvtk table2asn -i genbank_to_table/genbank_to_table.fsa -o temp.for.catching.errors \
    -s s_sourcefile.tsv -f genbank_to_table/genbank_to_table_fixed.tbl \
    -t template.sbt -c structured_comment.cmt | sed '1,3d' > suvtk_errors
rm temp.for.catching.errors.*

sed 's/\].*/]/' suvtk_errors | sort -u
grep "SEQ_FEAT.BadCDScomponentOverlapTRNA" suvtk_errors| grep -o '[ABC][0-9]*_[0-9]*_[0-9]*\.[0-9]*_[A-Z]\{2\}_[0-9]\{5\}_[a-z]\{3\}_[a-z]\{3\}_[dgh]' | sort -u > contigs_with_errors_trna
grep "SEQ_FEAT.NoStop" suvtk_errors| grep -o '[ABC][0-9]*_[0-9]*_[0-9]*\.[0-9]*_[A-Z]\{2\}_[0-9]\{5\}_[a-z]\{3\}_[a-z]\{3\}_[dgh]' | sort -u > contigs_with_errors_no_stop
grep "SEQ_FEAT.StartCodon" suvtk_errors| grep -o '[ABC][0-9]*_[0-9]*_[0-9]*\.[0-9]*_[A-Z]\{2\}_[0-9]\{5\}_[a-z]\{3\}_[a-z]\{3\}_[dgh]' | sort -u > contigs_with_errors_bad_start

#  tRNA entries are not allowed to (fully) overlap with CDS entries. (partial overlap is apparently ok)
cat genbank_to_table/genbank_to_table_fixed.tbl | \
    grep -v -P "67312\t67238\ttRNA" | \
    grep -v -P "67065\t66994\ttRNA" | \
    grep -v -P "23427\t23358\ttRNA" > genbank_to_table/genbank_to_table_fixed_twice.tbl

# CDSs that run off the edges need to have < before the left-column coordinate or > before the second-column coordinate (regardless of the strand, e.g. <999..1)
python $repo_location/scripts/HPC/Genbank_submission/fix_edge_genes.py genbank_to_table/genbank_to_table_fixed_twice.tbl contigs_with_errors_no_stop contigs_with_errors_bad_start genbank_to_table/genbank_to_table_fixed_thrice.tbl

# Now let's finally run this for good
suvtk table2asn -i genbank_to_table/genbank_to_table.fsa -o BPhage_genbank_submission \
    -s s_sourcefile.tsv -f genbank_to_table/genbank_to_table_fixed_thrice.tbl \
    -t template.sbt -c structured_comment.cmt 


# Second round...
conda activate viper_bphage
mkdir -p submission_2/genbank_to_table

seqkit grep -f submission_2/correct_these genbank_to_table/genbank_to_table.fsa  > submission_2/genbank_to_table/genbank_to_table2.fsa
head -n 1 s_sourcefile.tsv > submission_2/s_sourcefile2.tsv
grep -f submission_2/correct_these s_sourcefile.tsv >> submission_2/s_sourcefile2.tsv
python $repo_location/scripts/HPC/Genbank_submission/filter_feature_table.py genbank_to_table/genbank_to_table_fixed_thrice.tbl submission_2/correct_these submission_2/genbank_to_table/genbank_to_table_fixed_thrice2.tbl
head -n1 structured_comment.cmt > submission_2/structured_comment2.cmt
grep -f submission_2/correct_these structured_comment.cmt >> submission_2/structured_comment2.cmt

### MANUALLY CORRECT THE FEATURE TABLE submission_2/genbank_to_table/genbank_to_table_fixed_thrice2.tbl
