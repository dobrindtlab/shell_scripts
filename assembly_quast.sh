#!/bin/bash
set -e

### Script for automated assembly and QUAST quality control of forward and reverse fastq files.

if [ $# -lt 2 ]; then
	echo "Usage: $0 output_folder fastq-files-folder"
	exit 1
fi

[ -d "$1"assembly/ ] ||  mkdir "$1"assembly/;
[ -d "$1"quast_results/ ] ||  mkdir "$1"quast_results/;

for file in "$2"*_R1_001.fastq.gz;
do
prefix="${file%_R1_*}";
suffix="${prefix##*/}";
name="${suffix%_*}";
echo start assembly "$name";

### Assembly of raw reads with SPADES 3.5

spades.py --careful -t32 -k 21,31,55,77,99,127 -1 "$prefix"_R1_001.fastq.gz -2 "$prefix"_R2_001.fastq.gz -o "$1"assembly/"$name"/

[ -d "$1"assembly/all_fasta/ ] ||  mkdir "$1"assembly/all_fasta/;

### Exctract contigs larger than 1000bp

python long_seq.py "$1"assembly/"$name"/contigs.fasta "$1"assembly/all_fasta/"$name".fasta 1000;

### Analyse assembly quality with QUAST

echo start quast for "$name";
[ -d "$1"assembly/quast_results/ ] ||  mkdir "$1"assembly/quast_results/;
[ -d "$1"assembly/quast_results/"$name"/ ] ||  mkdir "$1"assembly/quast_results/"$name"/;

quast.py --threads 32 --est-ref-size 5500000 --labels "consensus, 21, 31, 55, 77, 99, 127" -o "$1"assembly/quast_results/"$name"/ "$1"assembly/"$name"/contigs.fasta "$1"assembly/"$name"/K21/final_contigs.fasta "$1"assembly/"$name"/K31/final_contigs.fasta "$1"assembly/"$name"/K55/final_contigs.fasta "$1"assembly/"$name"/K77/final_contigs.fasta "$1"assembly/"$name"/K99/final_contigs.fasta "$1"assembly/"$name"/K127/final_contigs.fasta;
done


### QUAST analysis of all consensus assemblies 
for file in "$1"assembly/all_fasta/; do
quast.py --threads 32 --min-contig 1000 --est-ref-size 5500000 -o "$1"quast_results/ "$1"assembly/all_fasta/*.fasta;
done

echo This Job is Done! AWESOME!!!

