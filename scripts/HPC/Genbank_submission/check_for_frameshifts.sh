
repo_location="$VSC_STAGING/BPhage"

cd $repo_location/scripts/HPC/Genbank_submission

python check_for_frameshifts.py ../genbank_to_table_fixed_thrice.tbl contig_list out_tbl

cat <(grep -A1 -B1 "contig" out_tbl) <(tail -n 1 out_tbl) | grep -v "contig" | grep -v "^--" | sort -u | \
    grep -v ">" | grep -v "<" > look_these_up_in_geneious