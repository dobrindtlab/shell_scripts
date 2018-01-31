#!/bin/bash
set -e

### Script for extracting a given DNA Sequence following the location and fasta_header names stored in the input file.


if [ $# -lt 3 ]; then
	echo "Usage: $0 input_fasta_file input_gene-location_file output_fasta_file"
	echo "The input_gene-location_file should be a tab-seperated file with 3 columns"
	echo "First column: faster_header, second colum: start location, third colum: stop location"
	exit 1
fi

while read fasta_header start_location stop_location; do
    echo "Write sequence for "$fasta_header""
    python /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/extract_DNA-seq_regions_from_fasta.py "$1" "$3" "$fasta_header" "$start_location" "$stop_location" 
 
done < "$2"
