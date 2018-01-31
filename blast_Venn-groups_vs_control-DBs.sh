#!/bin/bash
set -e

### script for the automatic analyse of extracted Venn group proteins of different Prokka runs.

if [ $# -lt 3 ]; then
	echo "Usage: $0 input_folder[extracted-proteins] inputfolder[BLAST DB] output-path[Blast]"
	exit 1
fi

[ -d "$3" ] ||  mkdir "$3"

for file in "$1"/*proteins.faa; do
	prefix="${file%*_proteins.faa}"
	name="${prefix##*/}"

### blast extracted proteins against control DBs.
	echo "---------- BLAST "$name"-specific proteins against Database ----------"
	blastp -task blastp -query "$file" -db "$2"/non-patho/non-patho -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$3""$name"_vs_non-patho.out
	blastp -task blastp -query "$file" -db "$2"/IPEC/IPEC -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$3""$name"_vs_IPEC.out
	blastp -task blastp -query "$file" -db "$2"/ExPEC/ExPEC -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$3""$name"_vs_ExPEC.out
	blastp -task blastp -query "$file" -db "$2"/IPEC_eae-pos/IPEC_eae-pos -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$3""$name"_vs_IPEC-eae-pos.out
	blastp -task blastp -query "$file" -db "$2"/IPEC_ST11_eae-pos/IPEC_ST11_eae-pos -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$3""$name"_vs_IPEC-ST11.out
	
done

echo "---------- done ----------"
