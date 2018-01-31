#!/bin/bash
set -e

### Skript to extract the DNA sequence of biomarker out of the Prokka ffn files.
### All Prokka results are stored in the ../Prokka/all_strains/output_Prokka/ folder!!

if [ $# -lt 3 ]; then
	echo "Usage: $0 input_file proteinortho_file output_folder"
	exit 1
fi
input="$1"
prefix="${input%*_IDs.txt}"
biomarker="${prefix##*/}"
proteinortho="$2"
output_folder="$3"

[ -d "$output_folder" ] ||  mkdir "$output_folder"
### 1. get the Protein ID out of the Biomarker_ID text file

while IFS='' read -r line || [[ -n "$line" ]]; do
    protein_ID="$line"
    echo "extract IDs for Biomarker "$protein_ID" and save them at "$output_folder""$protein_ID"_all-IDs.txt"

### 2. Take the Protein ID and search the proteinortho file for that ID. In the same line all other IDs are listed which are the same Protein.
### 3. Copy that line and extract the other Protein IDs	
	grep "$protein_ID" "$proteinortho" > "$output_folder""$protein_ID"_all-IDs.txt
	
### 3.1. need to delete all "*" and change all " " into "/n" and remove all empty lines!
	### change all "/t" into "/n" so that all information is sored in a new line:
	sed 's/\t/\n/g' "$output_folder""$protein_ID"_all-IDs.txt > "$output_folder""$protein_ID"_all-IDs_n.txt
	
	### remove all lines with only numbers
	awk '! /^[0-9]+$/' "$output_folder""$protein_ID"_all-IDs_n.txt > "$output_folder""$protein_ID"_all-IDs_no.txt
	### remove an occurence of 0.308 
	awk '! /^[0-9]+\.[0-9]+$/' "$output_folder""$protein_ID"_all-IDs_no.txt > "$output_folder""$protein_ID"_all-IDs_n2.txt
	### remove all lines with only a single character in it
	sed '/^.$/d' "$output_folder""$protein_ID"_all-IDs_n2.txt > "$output_folder""$protein_ID"_all-IDs.txt
	
	rm "$output_folder""$protein_ID"_all-IDs_n.txt
	rm "$output_folder""$protein_ID"_all-IDs_no.txt
	rm "$output_folder""$protein_ID"_all-IDs_n2.txt
	
### 4. Take each protein ID and extract the DNA sequence of each ID out of the all_strains/output_Prokka folder	

done < "$1"

for file in "$output_folder"*_all-IDs.txt; do
	pre="${file%_all-IDs.txt}"
	name="${pre##*/}"
	echo "I am saving all DNA seq for Biomarker "$name" right now"
	
	while IFS='' read -r line || [[ -n "$line" ]]; do
		Protein_ID="$line"
		strain="${Protein_ID%*_*}"
		echo "$line" > "$output_folder"samtools_input.txt
		xargs samtools faidx /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/all_250_strains/output_Prokka/"$strain"/"$strain".ffn < "$output_folder"samtools_input.txt >> "$output_folder""$biomarker"_"$name".fna
		rm "$output_folder"samtools_input.txt
	done < "$file"

done

