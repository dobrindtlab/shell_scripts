#!/bin/bash
set -e

### script for cross checking the Biomarker.
### 1. make BLAST DB for echt biomarker group
### 2. Blast each biomarker group against each BLAST DB. (cat all biomarker into one file??)
### 3. filter BLAST results as before. (same script? --> no adjust it!)


if [ $# -lt 1 ]; then
	echo "Usage: $0 input_folder"
	exit 1
fi


###  1. make BLAST DB for each biomarker group
[ -d "$1"BLAST_DBs/ ] ||  mkdir "$1"BLAST_DBs/

for file in "$1"*.faa; do
	prefix="${file%*_biomarker.faa}"
	name="${prefix##*/}"
	
	makeblastdb -in "$file" -dbtype 'prot' -title "$name" -out "$1"BLAST_DBs/"$name"
done

### 2. Blast each biomarker group against each BLAST DB. (--> cat all biomarker into one file!!)

cat "$1"*.faa > "$1"all_biomarker.fasta
[ -d "$1"BLAST_results/ ] ||  mkdir "$1"BLAST_results/

for file in "$1"*.faa; do
	prefix="${file%*_biomarker.faa}"
	name="${prefix##*/}"

### blast extracted proteins against control DBs.

	echo "---------- BLAST all biomarker against "$name" Database ----------"
	blastp -task blastp -query "$1"all_biomarker.fasta -db "$1"BLAST_DBs/"$name" -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$1"BLAST_results/all_biomarker_vs_"$name".out
done


### 3. filter BLAST results.
### extract all header of the all_biomarker.fasta file to get the index file for the R script

python extract_Fasta-header.py "$1"all_biomarker.fasta "$1"all_biomarker_index.txt
[ -d "$1"BLAST_results/filtered/ ] ||  mkdir "$1"BLAST_results/filtered/

for file in "$1"*ll_biomarker.faa; do
	prefix="${file%*_biomarker.faa}"
	name="${prefix##*/}"
	
	### take the name to choose the db as well as the index file! adjust the Rscript to take the index file into account!
	### quick and dirty without index file to mark the self hits!
	
	bash filtering_cross-check_blast-results.sh filtering_blast-results.R "$1"BLAST_results/ "$1"BLAST_results/filtered/ merge_biomarker_blast-results.R "$1"
#	"Usage: $0 path-to-R-script path-to-BLAST-results output-path[Blast] path-to-merge-R path-to-Protein-IDs"
	Rscript "$1" "$prefix"_sort.csv "$name" "$name"_80_80 "$db" "$db"_80/80 "$5""$Venn_group"_protein-IDs.txt "$3""$name"_best-hits.csv
done
