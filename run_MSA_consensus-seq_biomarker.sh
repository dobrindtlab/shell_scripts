#!/bin/bash
set -e

### script for multiple sequence alignment for biomarker primer design
### use clustalo for sequence alignment

for file in "$1"*_biomarker/*.fna; do
	prefix="${file%*.fna}"
	biomarker_name="${prefix##*/}"
	input="${prefix##*specific_biomarker/}"
	input_folder="${input%*/*}"
	
	echo "running script for biomarker: "$biomarker_name""
	[ -d "$1""$input_folder"/MSA/ ] ||  mkdir "$1""$input_folder"/MSA/
	clustalo -i "$file" --outfmt=clu -o "$1""$input_folder"/MSA/"$biomarker_name".aln

### use consambig out of the EMBOSS software package for detecting consenus sequence
	consambig -sequence "$1""$input_folder"/MSA/"$biomarker_name".aln -name "$biomarker_name" -outseq "$1""$input_folder"/MSA/"$biomarker_name".cons
done
