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
sinteractive --cpus-per-task=4 --mem=20g --gres=lscratch:50 
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
In the first step, we submit an R script (via Bash) that does the following
- Reads in a list of genes to check (one gene per line, header line = "gene_symbol")
- Reads in a csv of species to check from ensembl with four columns
	- Header row: genus_species,ensembl_name,common_name,group
	- genus_species ex. Xenopus_tropicalis
	- ensembl_name ex. xtropicalis (check biomart db list for this term, usually (always?) gspecies)
	- common_name ex. Tropical_clawed_frog
	- group ex. amphibian
- Uses biomaRt to search for each gene within each species database
- Finds the longest transcript for each gene found within each species
	- NOTE! If a fasta file with this name already exists, the script will skip querying ensembl for that species. This feature was added so that it is easier to add new species to the list later on and to make it easier to restart the program if it crashes because ensembl is temporarily down. This feature can be overriden using the --force flag during submission.
- Writes the peptide sequence to a user selected output folder
- Writes two swarm files, one to make the diamond indexes and one to run the diamond blastp searches
- Writes a script to submit those searches

```bash
sbatch 010_05_submit_diamond_prep.sh
```

### Step 010 > 020
This makes one diamond database per species. The source of the diamond database is the input file that was used for orthofinder (drawn from a copy stored in the 000_raw folder).
Next, the script runs a swarm of diamond blastp searches using the genes pulled from ensembl in the previous step as the query.
Each query is run against the appropriate species db. For example, the longest transcript for each requested gene from Xenopus tropicalis is queried against the list of X. tropicalis sequences used by Orthofinder.
Both of these swarm scripts are submitted via a bash script that includes a dependency so after all indexes are made, the blastp swarm gets automatically submitted
```bash
sbatch 020_05_run_diamond.sh
```
