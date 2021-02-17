## Libraries used
from argparse import ArgumentParser
import os
import re
import csv
import pathlib
from shutil import copyfile
from glob import glob
from collections import defaultdict, Counter


## Parser to handle user arguments
parser = ArgumentParser(description="Summarize orthogroups found across all genes.")

## Function to check if input directory exists
def input_dir_path(string):
		if os.path.isdir(string):
				return string
		else:
				raise NotADirectoryError(string)

## Input directory containg the gene-wise summary annotation files
parser.add_argument("-i", dest="in_dir", required=True,
										help="input directory containing annotation summaries from previous step (named 'by_gene/' by default)", metavar="Annotation_folder",
										type=input_dir_path)

## Input directory containg the orthogroup sequences from orthofinder
parser.add_argument("-d", dest="ofinder_dir", required=True,
										help="input directory orthogroup sequences from the orthofinder run", metavar="OG_seq_folder",
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

## Check if file location exists
def extant_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

## Location of csv file with species group info
parser.add_argument("-g", dest="group_info", required=True,
						help="Location of a CSV file containing the species names and the group they belong to (e.g., Xenopus_tropicalis,amphibian",
						metavar="Group_info.csv", type=extant_file)

## Parse arguments
args = parser.parse_args()


## Dictionary  connecting the species name to the species group
lib_to_group = {}
lib_to_species = {}

## Read in group info
with open(args.group_info) as f:
	reader = csv.reader(f)
	for line in reader:
		lib_to_species[line[2]]=line[0]
		lib_to_group[line[2]]=line[1]

print("Species groups processed.")



## Get all gene-wise csv summary files
gene_files = glob(args.in_dir + "*.csv")

OG_master_dict = {}

for gfile in gene_files:
	with open(gfile) as gf:
		## Get current gene name
		curr_gene = pathlib.Path(gfile).stem
		
		## Empty counter for OGs found for current gene
		OG_counter = Counter({'total':0})

		## Read gene csv
		reader = csv.reader(gf)

		## Skip header
		next(reader)


		for line in reader:

			## Fancy string to format with correct leading zeroes for CSV
			curr_OG = "OG" + str(int(line[-1])).zfill(7)

			## Count the number of times each unique OG shows up in a gene
			OG_counter[curr_OG] += 1
			OG_counter['total'] += 1

		## Store the counts
		OG_master_dict[curr_gene]=OG_counter

### SUMMARIZE ALL OGS FOUND ACROSS ALL GENES ###
## List of all OGs found for moving files later
all_ogs = []

out_file = args.out_dir + "/Orthogroup_agreement_summary.csv"

with open(out_file,'w') as f:
	writer = csv.writer(f)

	## Write headers
	headers = ['Gene','Orthogroup','N_agree','Pct_agree']
	writer.writerow(headers)

	## Write data for each OG found
	for gene_key in sorted(OG_master_dict.keys()):
		gene = gene_key
		OGs_found = OG_master_dict[gene].keys()

		for OG in OGs_found:
			if OG == "total":
				pass

			else:
				curr_OG_count = int(OG_master_dict[gene][OG])
				pct_agree = round(((float(curr_OG_count)/float(OG_master_dict[gene]['total']))*100),2)
				outline = [gene,OG,curr_OG_count,pct_agree]
				writer.writerow(outline)

				all_ogs.append("%s_%s.fa" % (gene,OG))

secondary_ogs = all_ogs.copy()

### SUMMARIZE TOP OG FOUND FOR EACH GENE ###
## List of top OGs found for moving files later
top_ogs = []

out_file = args.out_dir + "/Top_orthogroup.csv"

with open(out_file,'w') as f:
	writer = csv.writer(f)

	## Write headers
	headers = ['Gene','Orthogroup']
	writer.writerow(headers)

	## Write data for each OG found
	for gene_key in sorted(OG_master_dict.keys()):
		gene = gene_key
		OGs_dict = OG_master_dict[gene]

		## remove total from the dict
		OGs_dict.pop('total', None)

		## Then find the most common OG (must remove or total is always the top key)
		top_og = OGs_dict.most_common(1)[0][0]

		## Write top OG per gene
		outline = [gene,top_og]
		writer.writerow(outline)

		secondary_ogs.remove("%s_%s.fa" % (gene,top_og))
		top_ogs.append("%s_%s.fa" % (gene,top_og))




### COPY ORTHOGROUP SEQUENCES FILES
## Make landing folders for the orthogroup fasta files (if necessary)
top_og_dir = output_dir_path(args.out_dir + '/' + 'top_orthogroup_seqs')
secondary_og_dir = output_dir_path(args.out_dir + '/' + 'secondary_orthogroup_seqs')

## Copy Top OGs
for gene_og in top_ogs:
	og = gene_og.split('_')[-1]
	src_file = "%s/%s" % (args.ofinder_dir,og)
	dest_file = "%s/%s" % (top_og_dir,gene_og)
	copyfile(src_file, dest_file)

## Copy Supplmental/Secondary OGs
for gene_og in secondary_ogs:
	og = gene_og.split('_')[-1]
	src_file = "%s/%s" % (args.ofinder_dir,og)
	dest_file = "%s/%s" % (secondary_og_dir,gene_og)
	copyfile(src_file, dest_file)


### SUMMARIZE THE GROUP RICHNESS IN EACH OG
out_lines = []
for og in all_ogs:
	
	curr_out_line = []

	fileloc =  "%s/%s" % (top_og_dir,og)

	if os.path.exists(fileloc) == False:
		fileloc =  "%s/%s" % (secondary_og_dir,og)

	with open(fileloc) as f:
		gene = og.split('_')[0]
		og_id = og.split('_')[1]
		og_id = og_id.split('.')[0]
		total_seqs = 0
		libs_seen = []
		species_seen = []
		libs_seen_by_group = Counter({"Reptile":0, "Fish":0, "Amphioxus":0, "Cartilaginous_fish":0, "Tunicate":0, "Arthropod":0, "Hagfish":0, "Mollusc":0, "Bird":0, "Mammal":0, "Lamprey":0, "Amphibian":0, "Jellyfish":0})
		species_seen_by_group = Counter({"Reptile":0, "Fish":0, "Amphioxus":0, "Cartilaginous_fish":0, "Tunicate":0, "Arthropod":0, "Hagfish":0, "Mollusc":0, "Bird":0, "Mammal":0, "Lamprey":0, "Amphibian":0, "Jellyfish":0})
		
		for line in f:
			if line.startswith(">"):
				total_seqs += 1
				l = line.strip()
				lib_id = l.split("_")[-1]

				## Library counter
				if lib_id not in libs_seen:
					libs_seen.append(lib_id)

					group = lib_to_group[lib_id]
					libs_seen_by_group[group] += 1

				## Species counter
				sp = lib_to_species[lib_id]
				if sp not in species_seen:
					species_seen.append(sp)

					group = lib_to_group[lib_id]
					species_seen_by_group[group] += 1

	## Reformat for output
	curr_out_line.extend([gene, og_id, str(total_seqs), str(len(species_seen)),str(len(libs_seen))])

	for k in sorted(species_seen_by_group.keys()):
		curr_out_line.append(str(species_seen_by_group[k]))

	for k in sorted(libs_seen_by_group.keys()):
		curr_out_line.append(str(libs_seen_by_group[k]))

	out_lines.append(curr_out_line)

out_file = args.out_dir + "/Orthogroup_richness.csv"
with open(out_file,'w') as f:
	writer = csv.writer(f)

	headers = ['Gene','Orthogroup','Total_sequences','Unique_Species','Unique_Libraries']
	for k in sorted(species_seen_by_group.keys()):
		headers.append(k + "_Species")

	for k in sorted(libs_seen_by_group.keys()):
		headers.append(k + "_Libraries")

	writer.writerow(headers)
	writer.writerows(out_lines)

