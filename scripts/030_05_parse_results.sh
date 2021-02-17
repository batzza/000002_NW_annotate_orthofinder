#!/bin/bash
jobid1=$(python 030_10_identify_orthogroups.py -i ../data/020_diamond_output/ -o ../data/030_key_orthogroups/by_gene -s ../data/000_raw/SpeciesIDs.txt -q ../data/000_raw/SequenceIDs.txt -g ../data/000_raw/Group_info.csv -c ../data/000_raw/clusters_OrthoFinder_I1.5.txt_id_pairs.txt)
echo $jobid1
jobid2=$(sbatch --dependency afterany:$jobid1 python 030_15_summarize_orthogroups.py -i ../data/030_key_orthogroups/by_gene/ -d /data/VisionEvo/OrthoFinder/concat_unigene_108/OrthoFinder/Results_Feb04/Orthogroup_Sequences/ -o ../data/030_key_orthogroups/ -g ../data/000_raw/Group_info.csv)

