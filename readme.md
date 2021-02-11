# 000002_NW_annotate_orthogroups

**Codebase initiated**: 2021 February 11

**Author**: Zachary Batz

**Primary Investigator**: Noor White

## About

This project aims to annotate orthogroup output from a run of OrthoFinder.
The primary goal is to identify which gene families correspond to which orthogroup.
A secondary goal is to identify MAF-like genes in early chordates to check for funtional conservation of these TFs in a retinal context.

## Raw Data Source
Previously, we have run OrthoFinder to cluster unannotated originally found in /data/VisionEvo/OrthoFinder/Concat_unigene_pep/
 
OrthoFinder input was unigenes, concatenated by species with library-unique headers appended to individual sequence headers.
Concatenated species are named like this:
Mizuhopecten_yessoensis_C_denovo_unigenes.pep
 
Species for which we only had one sample are named like this:
Monodelphis_domestica_0_denovo_unigenes.pep 

All files were then renamed for OrthoFinder as Genus_species.fasta then run through OrthoFinder as detailed in /data/VisionEvo/OrthoFinder/how_to_run.txt

Because of the massive number of files produced by OrthoFinder, the data/000_raw folder for this project will not include the OrthoFinder output.
All OrthoFinder results will be accessed remotely from the original run directory: /data/VisionEvo/OrthoFinder/concat_unigene_108/OrthoFinder/Results_Feb04/

## Pipeline

### Setup Project Directory
```bash
python start_new_project.py NW annotate orthogroups
```

### Initiate Conda Env
```bash
conda create -n 000002_NW
conda activate 000002_NW
```

### Install Necessary Modules
TBD

### Setup Git Repo
```bash
git init
git add scripts/
git add readme.md
git commit -m "Subject: Abc" -m "Body: xyz"
git remote add origin https://github.com/batzza/000002_NW_annotate_orthofinder.git
git branch -M main
git push -u origin main
```

### Step 000 > 010
First we will need to get the appropriate reference sequences, maybe from RefSeq or from Ensembl?
