#!/bin/bash
set -e

### bash script to call the R analysis and alter the BLAST output before!


if [ $# -lt 5 ]; then
	echo "Usage: $0 path-to-R-script path-to-BLAST-results output-path[Blast] path-to-merge-R path-to-Protein-IDs"
	exit 1
fi

[ -d "$3" ] ||  mkdir "$3"

for file in "$2"*.out; do
	prefix="${file%*.out}"
	name="${prefix##*/}"
	db="${name##*_}"
	Venn_group="${name%*_vs*}"
	
	### change from tab delimited file to ";" delimited file:
	sed 's/\t/;/g' "$file" > "$prefix".csv
	
	### delete all line with "#" at the start --> not important information and could not read by R
	sed '/^#/ d' "$prefix".csv > "$prefix"_sort.csv
	
	### input prepared file into R for analysis:
	### usage: R R_schript BLAST_results.csv Name_der_BLAST_kombi output-PATH
	
	if [[ -s "$prefix"_sort.csv ]] ### checks if file is empty or not
		then
			echo "run Rscript for "$Venn_group" with database "$db"."	
			Rscript "$1" "$prefix"_sort.csv "$name" "$name"_80_80 "$db" "$db"_80/80 "$5""$Venn_group"_protein-IDs.txt "$3""$name"_best-hits.csv "$3""$name"_biomarker-pair.txt
			echo "done"
		else
			echo "file "$prefix"_sort.csv is empty. Need to work around"
			echo "Protein_ID,"$db","$db"_80/80" > "$3""$name"_best-hits.csv 
	fi		
			
	#### if all three BLAST scores for one VENN Group are completed, merge them!
	if [ -e "$3""$Venn_group"_vs_non-patho_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_IPEC_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ExPEC_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_IPEC-ST11_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_IPEC-eae-pos_best-hits.csv ]
		then 
			echo "------------all BLAST best hits done------------"
			echo "------------merge BLAST best hits---------------"
			### use a small Rscript to merge the five .csv files for each database together!
			Rscript "$4" "$3""$Venn_group"_vs_non-patho_best-hits.csv "$3""$Venn_group"_vs_IPEC_best-hits.csv "$3""$Venn_group"_vs_ExPEC_best-hits.csv  "$3""$Venn_group"_vs_IPEC-eae-pos_best-hits.csv "$3""$Venn_group"_vs_IPEC-ST11_best-hits.csv "$3""$Venn_group"_merged_best-hits.csv
			echo "------------"$Venn_group" BLAST results are merged-----------"
		else  
			echo "------------at least one BLAST best hit in "$Venn_group" still missing------------"
	fi
	### remove file to save disk space and keep clarity
	rm "$prefix".csv
	rm "$prefix"_sort.csv
done



	
