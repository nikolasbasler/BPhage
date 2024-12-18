##############################################################################
###### INCREASE CONFIDENCE THAT OUR PHAGES ARE NOT BACTERIAL CONTAMINATION
##############################################################################

# Jelle wanted to see if the cluster members are all of similar length, to have more confidence that they are not bacterial contaminants.
# I don't think this is a good way of testing this, because the clusters will contain contigs of various length anyway.

cd $VSC_STAGING/BPhage/output/
grep -f ../data/core_contigs.txt <(zcat bphage_ALL_1kb_cross_95-85_clusters.tsv.gz) | cut -f2 > temp.core_clusters

rm -f temp.cluster_histograms
while read line; do
    echo $line | cut -d"," -f1 >> temp.cluster_histograms
    echo $line | tr ',' '\n' | cut -d "_" -f4 | $VSC_STAGING/histogram.awk -v bw=5000 >> temp.cluster_histograms
    echo "" >> temp.cluster_histograms
done < temp.core_clusters
# -> The cluster representatives are always the longest by at least some margin. Most members are always below 5k.

# Instead, check if genomad flagged any host regions that we cut off during per-sample clustering (indicated by the V or H in the length string)

zcat bphage_ALL_1kb_cross_95-85_clusters.tsv.gz | cut -f2 | tr ',' '\n' | grep -e "[0-9]H" -e "V" \
    > temp.cluster_members_with_host_regions
grep -f temp.cluster_members_with_host_regions <(zcat bphage_ALL_1kb_cross_95-85_clusters.tsv.gz) | cut -f1 \
    > temp.reps_with_members_with_host_regions
wc -l temp.reps_with_members_with_host_regions

cat bphage_ALL_1kb_phages.csv bphage_ALL_1kb_picobirna.csv bphage_ALL_1kb_unclassified_viruses.csv | \
    grep -f temp.reps_with_members_with_host_regions | cut -d"," -f1 \
    > temp.phage_clusters_with_members_with_host_regions
wc -l temp.phage_clusters_with_members_with_host_regions
grep -e "[0-9]H" -e "V" temp.phage_clusters_with_members_with_host_regions | wc -l

grep -f ../data/core_contigs.txt temp.reps_with_members_with_host_regions > temp.core_clusters_with_members_with_host_regions
wc -l temp.core_clusters_with_members_with_host_regions
grep -e "[0-9]H" -e "V" temp.core_clusters_with_members_with_host_regions | wc -l

# -> Only 44 of the 1,216,077 total contigs (members and reps) were flagged (35 of the 287,278 representatives),
#       6 of them are part of the 2,346 contigs classified as phages (4 of them are representatives),
#       1 of the 97 core phages (not a representative)

zcat bphage_ALL_1kb_cross_95-85_clusters.tsv.gz | grep -f temp.core_clusters_with_members_with_host_regions | cut -f2 | tr ',' '\n' \
    > temp.core_contig_with_host_cluter_members
wc -l temp.core_contig_with_host_cluter_members
grep -e "[0-9]H" -e "V" temp.core_contig_with_host_cluter_members > temp.the_member_with_host
wc -l temp.the_member_with_host

grep -A5 -f temp.core_clusters_with_members_with_host_regions temp.cluster_histograms
cat temp.the_member_with_host

# -> The 1 core phage that has members with host regions is from a large cluster of 577 contigs. Only one of them is flagged.
#   It was stripped from 13.6k down to 5.5k. Looking at the classification table in R, the representative has no host genes,
#   but is 50.88% complete and happens to be the one that occurs in all hives. I still don't think it needs any intervention.

##############################################################################
###### 
##############################################################################
