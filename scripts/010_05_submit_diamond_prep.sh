#!/bin/bash
module load R/4.0
Rscript 010_10_prep_diamond.r --gene_list ../data/000_raw/gene_list.txt --species_list ../data/000_raw/species_list.txt -o ../data/010_diamond_input --hit_dir ../data/020_diamond_output --threads 2
