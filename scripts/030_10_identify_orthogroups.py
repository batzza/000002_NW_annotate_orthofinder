## Libraries used
from argparse import ArgumentParser
import os
import re
import csv
import pathlib
from glob import glob
from collections import defaultdict


## Parser to handle user arguments
parser = ArgumentParser(description="Identify orthogroups of interest by parsing results from the orthofinder annotation script.")

## Input directory argument and functino function to check if input directory exists
def input_dir_path(string):
		if os.path.isdir(string):
				return string
		else:
				raise NotADirectoryError(string)

parser.add_argument("-i", dest="in_dir", required=True,
										help="input directory containing output of diamond blastp searches in default BLAST tabular format", metavar="Diamond_outptu",
										type=input_dir_path)

## Output directory argument and function to make it if it does not already exist
def output_dir_path(string):
		if os.path.isdir(string):
				return string
		else:
				pathlib.Path(string).mkdir(parents=True, exist_ok=True)
				return string

parser.add_argument("-o", dest="out_dir", required=True,
										help="Directory for output files", metavar="OUTPUT_DIRECTORY",
										type=output_dir_path)

## Location of SpeciesIDs.txt file made by OrthoFinder
## Check if file location exists
def extant_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

parser.add_argument("-s", dest="speciesIDs", required=True,
						help="Location of SpeciesIDs.txt file made by OrthoFinder",
						metavar="SpeciesIDs.txt", type=extant_file)

## Location of SequenceIDs.txt file made by OrthoFinder
parser.add_argument("-q", dest="sequenceIDs", required=True,
						help="Location of SequenceIDs.txt file made by OrthoFinder",
						metavar="SequenceIDs.txt", type=extant_file)

## Location of csv file with species group info
parser.add_argument("-g", dest="group_info", required=True,
						help="Location of a CSV file containing the species names and the group they belong to (e.g., Xenopus_tropicalis,amphibian",
						metavar="Group_info.csv", type=extant_file)

## Location of csv file with species group info
parser.add_argument("-c", dest="orthogroups", required=True,
						help="Location of a clusters info file from OrthoFinger output",
						metavar="clusters_OrthoFinder_I1.5.txt_id_pairs.txt", type=extant_file)

## Parse arguments
args = parser.parse_args()

## Dictionary  connecting the species name to the species group
## Key is species name (Genus)_(species); value is number in OF
species_group_dict = {}

## Read in group info
with open(args.group_info) as f:
	reader = csv.reader(f)
	for line in reader:
		species_group_dict[line[0]]=line[1]

print("Species groups processed.")

## Dictionary  connecting the species name to the species ID
## Key is species name (Genus)_(species); value is number in OF
species_id_dict = {}

## Read in species id info
with open(args.speciesIDs) as f:
	for line in f:
		spID = line.split(':')[0]
		spName = line.split(':')[1]
		spName = spName.split('.')[0]
		spName = spName[1:]
		species_id_dict[spName] = spID

print("Species IDs processed.")

## Dictioanry connecting trinity sequence name (orthofinder input) with orthofinder alias
## Key is trinity seq name; value is orthofinder alias seq name
sequence_id_dict = {}
with open(args.sequenceIDs) as f:
	for line in f:
		seqID = line.split(':')[0]
		seqName = line.split(':')[1]
		seqName = seqName[1:]
		sequence_id_dict[seqName.strip()] = seqID

print("Sequence IDs processed.")


## Dictionary connecting OF sequence ID aliases to orthogroups
## Key is sequence ID alias, value is orthogroup
og_dict = {}
with open(args.orthogroups) as f:
	n=0
	ignoreFlag = 1
	for l in f:
		if ignoreFlag == 1:
			if l.startswith("begin"):
				ignoreFlag = 0

		else:
			n+=1
			if n % 500000 == 0:
				print("%s clusters lines read..." % n)
			lclean = l.strip()
			line = lclean.split(" ")
			

			if l[0].isdigit():
				curr_og = line[0]
				hits = line[1:]
				
			else:
				hits = line

			for seq in hits:
				if seq != "$" and len(seq)>1:
					og_dict[seq]=curr_og


print("Orthogroup members processed.")

## Find all the diamond top hit files
hit_files = glob(args.in_dir+"*_top_hits.txt")

## Dictionary where the key is the gene and the results are a list of top hits for that gene (including the alignment stats)
alignment_dict = defaultdict(list)

for h in hit_files:
	print("Processing %s." % h)
	curr_sp = os.path.basename(h)
	curr_sp = re.sub("_top_hits.txt","",curr_sp)
	
	try:
		curr_grp = species_group_dict[curr_sp]
	except KeyError:
		curr_grp = "NA"

	with open(h) as f:
		for l in f:

			lclean = l.strip()
			line = lclean.split('\t')
			gene = line[0]
			hit_deets = [curr_sp]
			hit_deets.append(curr_grp)

			hit_deets.append(line[1])
			hit_deets.append(line[2])
			hit_deets.append(line[3])
			hit_deets.append(line[10])
			hit_deets.append(line[11])

			of_id = sequence_id_dict[line[1]]
			hit_deets.append(of_id)

			og_id = og_dict[of_id]
			hit_deets.append(og_id)

			## Feb 16 2021 -- there is something a little strange about the results
			## search results from Diamond are duplicated 2-3x times in most _top_hits.txt files
			## The results are exactly identical so I think somehow the files just got appended instead of overwritten at some point?
			## The end result is that the same results get printed multiple times to the output CSV files
			## Not the end of the world but a little annoying
			## Rather than re-run the diamond searches, I am adding a test to check if the result list already exists in the dictionary entry for the current gene

			## Check if the list already exists
			## If the key exists and the current hit detail list is already a stored value, pass
			## If the key exists but the current hit detail list is NOT already stored, store it
			## If the key does not exists, store the key and the value
			try:
				if hit_deets in alignment_dict[gene]: 
					pass
				else: 
					alignment_dict[gene].append(hit_deets)

			except KeyError:
				alignment_dict[gene].append(hit_deets)


print("Writing results file.")
print(alignment_dict.keys())
for k in alignment_dict.keys():
	outfile = args.out_dir + '/' + k + ".csv"
	print(outfile)
	with open(outfile,'w') as f:
		result_list = alignment_dict[k]

		writer = csv.writer(f)
		writer.writerow(['Ref_species','Ref_group','Trinity_id','Aln_pct','Aln_len','Eval','Bitscore','OrthoFinder_seq','Orthogroup'])
		writer.writerows(result_list) 

print("Done.")

