#!/bin/bash
set -e

### Script for automated BLASTp anaylsis

if [ $# -lt 4 ]; then
	echo "Usage: $0 path_to_faa-files name-of-DB proteins-to-BLAST.faa output-path"
	exit 1
fi

### create database of all protein sequences:
cat "$1"*.faa > "$1"all_seq.fasta

### create BLASTp database for all sequences:
makeblastdb -in "$1"all_seq.fasta -dbtype 'prot' -title "$2" -out "$1"blast-db/"$2"

### blast proteins against BLASTp database.
blastp -task blastp -query "$3" -db "$1"blast-db/"$2" -evalue 1e-5 -num_threads 32 -parse_deflines -perc_identity 90 -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$4""$2"_90.out



