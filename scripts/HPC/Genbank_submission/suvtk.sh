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
python3 $repo_location/scripts/HPC/Genbank_submission/filter_gbk_and_remove_locus_tag.py s_all_genomes.gbk s_genomes.list s_genomes.gbk

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

# Add the missing info for repeat_region to the table (it's in pharokka's minced gff output, unfortunately under different names...)
sed 's/48340\t48528\trepeat_region/48340\t48528\trepeat_region\n\t\t\trpt_type\tdirect\n\t\t\trpt_family\tCRISPR\n\t\t\trpt_unit_seq\tAAAGGTGTCTGGGAATATTCAAATAGCATTG/g' genbank_to_table/genbank_to_table_fixed_thrice.tbl > temp.fixed_repeat_regions
sed -i 's/28834\t29044\trepeat_region/28834\t29044\trepeat_region\n\t\t\trpt_type\tdirect\n\t\t\trpt_family\tCRISPR\n\t\t\trpt_unit_seq\tTTTCTAAACCGCTTATGCAGCGGTGAACAGT/g' temp.fixed_repeat_regions
sed -i 's/32917\t33244\trepeat_region/32917\t33244\trepeat_region\n\t\t\trpt_type\tdirect\n\t\t\trpt_family\tCRISPR\n\t\t\trpt_unit_seq\tTTTCTAAACCGCCTATTCGGTGGTAAAC/g' temp.fixed_repeat_regions
sed -i 's/9804\t10076\trepeat_region/9804\t10076\trepeat_region\n\t\t\trpt_type\tdirect\n\t\t\trpt_family\tCRISPR\n\t\t\trpt_unit_seq\tATGTTCCCTGTATGCACAGGGATAAACCG/g' temp.fixed_repeat_regions

# Fix the cases where the CDS should be partial but isn't for some reason:
sed -i 's/1\t1404\tCDS/<1\t1404\tCDS/g' temp.fixed_repeat_regions
sed -i 's/3\t1331\tCDS/<3\t1331\tCDS/g' temp.fixed_repeat_regions 
sed -i 's/^2\t334\tCDS/<2\t334\tCDS/g' temp.fixed_repeat_regions 
sed -i 's/^3964\t3731\tCDS/<3964\t3731\tCDS/g' temp.fixed_repeat_regions
sed -i 's/^4329\t3772\tCDS/<4329\t3772\tCDS/g' temp.fixed_repeat_regions

# Fixing clipped CDSs

python $repo_location/scripts/HPC/Genbank_submission/fix_coordinates_of_edge_genes.py temp.fixed_repeat_regions temp.fixed_clipped_cds

# Adding tRNA AA (I used a different scheme to shorten the names for pharokka...)
sed 1d genbank_investigation/second_round/missing_tRNA_AA | awk '{print $2}' | cut -d "|" -f2 | sed 's/:/\t/g' | sed 's/-/\t/g' > temp.trna.coords
cut -f1 temp.trna.coords | awk -F "_" '{print "NODE_"$1"_"$4"_"$5"_"$6"_"$7"_"$8}' > temp.trna.pharokka_names
paste temp.trna.coords temp.trna.pharokka_names | sed 's/c\([0-9].*\)/\1/g' > temp.trna.names_and_coords
while read line; do
    name=$(echo "$line" | awk '{print $4}')
    real_name=$(echo "$line" | awk '{print $1}')
    coord_one=$(echo "$line" | awk '{print $2}')
    coord_two=$(echo "$line" | awk '{print $3}')
    if [ $coord_one -gt $coord_two ]; then
        scan_coord_one=$coord_two
        scan_coord_two=$coord_one
    else
        scan_coord_one=$coord_one
        scan_coord_two=$coord_two
    fi
    raw_line=$(awk -v name="$name" -v coord_one="$scan_coord_one" -v coord_two="$scan_coord_two" '$1==name && $3!="exon" && $4==coord_one && $5==coord_two {print $9}' $repo_location/output/annotation/pharokka_bphage_and_others/trnascan_out.gff)
    anticodon=$(echo $raw_line | cut -d";" -f2 | rev | cut -d "-" -f1 | cut -c -3 | rev)
    aminoacid=$(echo $raw_line | cut -d";" -f2 | rev | cut -d "-" -f1 | cut -c 4- | rev)   
    echo -e "$real_name\t$coord_one\t$coord_two\t$anticodon\t$aminoacid"
done < temp.trna.names_and_coords > trna_anticodons_and_aminoacids

python $repo_location/scripts/HPC/Genbank_submission/fix_tRNAs.py temp.fixed_clipped_cds trna_anticodons_and_aminoacids genbank_to_table/genbank_to_table_fixed_thrice_repeats_clips_trnas.tbl


# Now let's finally run this for good
# suvtk table2asn -i genbank_to_table/genbank_to_table.fsa -o BPhage_genbank_submission \
#     -s s_sourcefile.tsv -f genbank_to_table/genbank_to_table_fixed_thrice.tbl \
#     -t template.sbt -c structured_comment.cmt 

suvtk table2asn -i genbank_to_table/genbank_to_table.fsa -o BPhage_genbank_submission_updated \
    -s s_sourcefile.tsv -f genbank_to_table/genbank_to_table_fixed_thrice_repeats_clips_trnas.tbl \
    -t template.sbt -c structured_comment.cmt 

rm temp.*
