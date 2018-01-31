#!/bin/bash
set -e

### bash script to identify ORFs and analyse distribution with Proteinortho5 and display the results groupwise in venn diagrams.
### Requirements for usage:
###		1. 	copy fasta files in new folder fasta_files/ .
###		2. 	create group_file.tsv in input_folder with all sequencenames group in max 4 columns (first colum (header) is name of the group).
### 		IMPORTANT: within group_file.tsv all genome names must have a ".faa" ending!
###		3. 	Usage: bash Prokka_ORF_finding.sh "Projekt" "input_folder"
###			"input_folder" must contain the hard path! (not relative!)	
### 	4.	run script in same folder then accompanying perl scripts (https://github.com/aleimba/bac-genomics-scripts).


if [ $# -lt 2 ]; then
	echo "Usage: $0 Project input_folder[absolut path]"
	exit 1
fi

project="$1"
input_folder="$2"

[ -d "$input_folder"output_Prokka/ ] ||  mkdir "$input_folder"output_Prokka/

shopt -s nullglob
for file in "$input_folder"fasta_files/*.fasta; do
	prefix="${file%*.fasta}"
	mv "$file" "$prefix".fna
done
	
for file in "$input_folder"fasta_files/*.fna; do
	prefix="${file%*.fna}"
	name="${prefix##*/}"

	perl prokka --outdir "$input_folder"output_Prokka/"$name"/ --force --prefix "$name" --addgenes --locustag "$name" --centre MS --compliant --genus Escherichia --species coli --strain "$name" --usegenus --cpus 30 "$prefix".fna
	
	[ -d "$input_folder"output_Prokka/"$name"/origin_files/ ] ||  mkdir "$input_folder"output_Prokka/"$name"/origin_files/
	cp "$input_folder"output_Prokka/"$name"/"$name".faa "$input_folder"output_Prokka/"$name"/origin_files/"$name"_origin.faa
	mv "$input_folder"output_Prokka/"$name"/"$name".gff "$input_folder"output_Prokka/"$name"/origin_files/"$name"_origin.gff
	
	perl cds_extractor.pl -i "$input_folder"output_Prokka/"$name"/"$name".gbk -p
	
	cp "$input_folder"output_Prokka/"$name"/"$name".faa "$input_folder"output_Prokka/"$name".faa 
	
	perl bp_genbank2gff3.pl -r "$input_folder"output_Prokka/"$name"/"$name".gbk -o "$input_folder"output_Prokka/"$name"/
	
	mv "$input_folder"output_Prokka/"$name"/"$name".gbk.gff  "$input_folder"output_Prokka/"$name"/"$name".gff
	
	for gff in "$input_folder"output_Prokka/"$name"/*.gff; do
		cp "$gff" "$input_folder"output_Prokka/"$name"/origin_files/"$name"_with-fasta.gff
		perl -ne 'last if (/\#\#FASTA/);print;' "$gff" > "$gff"-no-fasta
		mv "$gff"-no-fasta "$gff"
	done
done

[ -d "$input_folder"Protein_ortho/ ] || mkdir "$input_folder"Protein_ortho/
cd "$input_folder"Protein_ortho/

proteinortho5.pl -selfblast -graph -singles -synteny -verbose --keep -identity=70 -cov=70 -blastParameters='-use_sw_tback' -project="$project"_i70_c70 -clean -cpus=32 "$input_folder"output_Prokka/*/*.faa

cp "$project"*.proteinortho ../

[ -d "$input_folder"ORF_results/ ] || mkdir "$input_folder"ORF_results/

for proteinmatrix in "$input_folder"*.proteinortho; do
	protprefix="${proteinmatrix%*.proteinortho}"
	protname="${protprefix##*/}"
	perl_script="po2group_stats.pl"
	group_file=""$input_folder"group_file.tsv"
	ortho_input="$proteinmatrix"
	echo Orthologmatrix wird verglichen für Projekt: "$protname"

	[ -d "$input_folder"ORF_results/"$protname"/ ] || mkdir "$input_folder"ORF_results/"$protname"/
	perl "$perl_script" -i "$ortho_input" -d "$input_folder"output_Prokka/ -g "$group_file" -r "$input_folder"ORF_results/"$protname"/strict_100/ -cut_i 1 -cut_e 0.0 -b -p -s -u > "$input_folder"ORF_results/"$protname"/strict_"$protname"_overall_stats.tsv
done
	


	
