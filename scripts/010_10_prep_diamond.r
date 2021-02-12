library(biomaRt)
suppressPackageStartupMessages(library(dplyr))
library(getopt)

## Describe the expected command line arguments
spec <- matrix(c(
  # long name  short name  mask  type          description(optional)
  # ---------  ----------  ----  ------------  ---------------------
  'gene_list','g', 1, 'character','text file containing the genes of interest, one gene per line',
  'species_list','s',1,'character','text file containing the species of interest, one species per line ',
  'out_dir','o',1,'character','directory where fastas with ensembl genes will go',
  'hit_dir','d',1,'character','directory where the diamond hits will go',
  'threads','t',0,'integer','number of threads to request for diamond swarms (default = 2)',
  'force','f',0,  'logical','By default, this program will skip creating any fasta files that already exist. Use -f to force them to be regenerated',
  'help','h',0,  'logical','show this help message'
), byrow=TRUE, ncol=5);

## parse the command line
opt <- getopt(spec);

## show help if requested
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE));
  q();
}

## set defaults
if ( is.null(opt$threads) )    { opt$threads = 2 }

## Make output directories
outdir_main <- opt$out_dir
outdir_db <- paste0(outdir_main,"/dbs/")
outdir_qs <- paste0(outdir_main,"/queries/")
dir.create(file.path(outdir_main), showWarnings = FALSE)
dir.create(file.path(outdir_db), showWarnings = FALSE)
dir.create(file.path(outdir_qs), showWarnings = FALSE)
dir.create(file.path(opt$hit_dir), showWarnings = FALSE)

## Load the ensembl mart
ensembl <- useEnsembl(biomart = "ensembl")

## Store current version ID
version <- tail(strsplit(listMarts()[1,"version"]," ")[[1]],1)

## Read in the list of species to work with
sp_list <- read.csv(opt$species_list)

## Read in the list of genes to pull
gene_list <- read.csv(opt$gene_list)

## List of something... tbd
anno_list <- c()


for (sp in sp_list$ensembl_name){
  ## Check if file exists already
  fname <- sprintf("/v%s_%s_longest_gene.fasta", version, sp)
  outname <- paste0(outdir_qs,fname)
  
  ## Only attempt to query ensembl if the fasta file does not exist or the user has included the --force option
  if ((!file.exists(outname)) || ((!file.exists(outname)) && opt$force == TRUE)){
    
  ## Prepare ensembl to query the current species
  dset <- paste0(sp, "_gene_ensembl")
  ensembl <- useDataset(dataset = dset, mart = ensembl)
  
  ## Get Ensembl sequences -- all transcripts for genes in list for current species
  ensembl_hits<-getSequence(id=gene_list$gene_symbol,
                            type="external_gene_name",
                            seqType = "peptide",
                            mart=ensembl)

  ## Remove sequences that have no associated transcript mRNA  
  ensembl_hits<-ensembl_hits[which(ensembl_hits$peptide != "Sequence unavailable"),]
  
  ## Sort by gene name and make all gene names upper case
  ensembl_hits<- ensembl_hits %>% arrange(external_gene_name)
  ensembl_hits$external_gene_name <- toupper(ensembl_hits$external_gene_name)
  
  ## Keep just the longest transcript for each gene
  ensembl_hits <- ensembl_hits %>% group_by(external_gene_name) %>% slice_max(length(peptide), n = 1) 
  
  ## dplyr adds some attributes and modifies the character vectors in a way that the output exportFasta balks at
  ## Messed around a little and couldn't figure out how to fix this in place so my hack solution is to
  ## just move the two columns I need to a new dataframe and export that -- works fine that way
  out_seqs <- data.frame(peptide = ensembl_hits$peptide)
  out_seqs$external_gene_name <- ensembl_hits$external_gene_name

  fname <- sprintf("/v%s_%s_longest_gene.fasta", version, sp)
  outname <- paste0(outdir_qs,fname)
  exportFASTA(out_seqs,outname)
  }
  
  anno_list <- c(anno_list, fname)
}



anno_tbl <- tibble(genus_species = sp_list$genus_species, anno_list)
anno_tbl$index_db <- sprintf("diamond makedb --in ../data/000_raw/OrthoFinder_input/%s.fasta --db %s%s.db --threads %i", anno_tbl$genus_species, outdir_db, anno_tbl$genus_species, opt$threads)
anno_tbl$diamond_run <- sprintf("diamond blastp --db %s%s.db --query %s%s --mid-sensitive --max-target-seqs 1 --unal 0 --threads %i --out %s/%s_top_hits.txt", outdir_db, anno_tbl$genus_species, outdir_qs, anno_tbl$anno_list, opt$threads, opt$hit_dir, anno_tbl$genus_species)

## Write submission bash script
outname_submit <- "020_05_run_diamond.sh"

## Write swarm file to make indexes for Diamond
outname_db <- "020_10_make_indexes.swarm"
write.table(anno_tbl$index_db, outname_db, row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)

## Write swarm of jobs to run Diamond search
outname_search <- "020_15_run_diamond.swarm"
write.table(anno_tbl$diamond_run, outname_search, row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)




## This section writes a swarm submission where the second job is automatically submitted after the first job is completed
runLines <- c()
runLines <- c(runLines, sprintf("#!/bin/bash"))
runLines <- c(runLines, sprintf("jobid1=$(swarm -f %s --module diamond)", outname_db))
runLines <- c(runLines, sprintf("echo $jobid1"))
runLines <- c(runLines, sprintf("jobid2=$(swarm -f %s --module diamond --dependency afterany:$jobid1)", outname_search))

fileConn<-file(outname_submit)
writeLines(runLines, fileConn)
close(fileConn)
