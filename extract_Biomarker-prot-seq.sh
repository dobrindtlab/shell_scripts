#!/bin/bash
set -e


### Skript to extract protein sequences from putative biomarker index

if [ $# -lt 1 ]; then
	echo "Usage: $0 input_folder(Prokka-analysis)"
	exit 1
fi

for file in "$1"BLAST/results/biomarker/*.txt; do
	prefix="${file%*_biomarker.txt}"
	name="${prefix##*/}"
	
	echo "------ extracting biomarker for "$name" -------"
	xargs samtools faidx "$1"/extracted-proteins_strict_80-40/"$name"_proteins.faa < "$file" > "$1"BLAST/results/biomarker/"$name"_biomarker.faa
	
done
