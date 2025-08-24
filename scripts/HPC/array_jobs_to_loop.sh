#!/bin/bash -l

repo_location="$VSC_STAGING/BPhage"

cd $repo_location/scripts/HPC/

grep -l "SLURM_ARRAY_TASK_ID" *slrm | grep -v "slurm_template.slrm" > temp.array_scripts

while read script; do
    loop_file=$(grep "SLURM_ARRAY_TASK_ID" $script | grep -v "^#" | cut -d"/" -f3 | cut -d" " -f1)
    rep_loc_line=$(grep "repo_location=" $script)

    echo $rep_loc_line > temp.script
    echo "while read line; do" >> temp.script
    grep -v "SLURM_ARRAY_TASK_ID" $script | grep -v "repo_location=" >> temp.script
    echo "done < \$repo_location/data/$loop_file" >> temp.script

    mv temp.script $script

done < temp.array_scripts

sed -i 's/done.*/done <(cut -f1 $repo_location\/data\/other_datasets_SRA_accessions.tsv)/g' additional_datasets_mapping_with_unpaired.slrm 

# Clean up
rm temp.array_scripts
