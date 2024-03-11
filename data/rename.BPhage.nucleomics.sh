#! /bin/bash

# NOTE: Sample 4295 failed in sequencing run R4317 and was re-run 
# in R4341.

# The sequencing was done using the tubes ID as library names. This
# script renames the files according to the sample ID, which is more
# meaningful.

while read -r line;
do
	nuc=$(echo $line | awk '{print $1}')
	me=$(echo $line | awk '{print $2}')

	mv ${nuc}*R1*fastq.gz ${me}.R1.fastq.gz
        mv ${nuc}*R2*fastq.gz ${me}.R2.fastq.gz

done < rename.BPhage.nucleomics.scheme
