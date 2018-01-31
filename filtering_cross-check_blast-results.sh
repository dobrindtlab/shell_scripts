#!/bin/bash
set -e

### bash script to call the R analysis and alter the BLAST output before!
### dependencies: R Scripts: - filtering_cross-check-blast-results.R
###							 - merge_biomarker_blast-results.R


if [ $# -lt 5 ]; then
	echo "Usage: $0 path-to-filtering-R-script path-to-BLAST-results output-path[Blast] path-to-merge-R-script path-to-Protein-IDs"
	exit 1
fi

[ -d "$3" ] ||  mkdir "$3"

for file in "$2"*.out; do
	prefix="${file%*.out}"
	name="${prefix##*/}"
	db="${name##*vs_}"
	Venn_group="${name%*_vs*}"
	
	### change from tab delimited file to ";" delimited file:
	sed 's/\t/;/g' "$file" > "$prefix".csv
	
	### delete all line with "#" at the start --> not important information and could not read by R
	sed '/^#/ d' "$prefix".csv > "$prefix"_sort.csv
	
	### input prepared file into R for analysis:
	### usage: R R_schript BLAST_results.csv name_of_BLAST_combination output-PATH
	### Problem: If no hits were found in the VENN analysis, files are empty! Problem with R!
	
	if [[ -s "$prefix"_sort.csv ]] ### checks if file is empty or not
		then
			echo "run Rscript for "$Venn_group" with database "$db"."	
			Rscript "$1" "$prefix"_sort.csv "$name" "$name"_80_80 "$db" "$db"_80/80 "$5""$Venn_group"_index.txt "$3""$name"_best-hits.csv "$3""$name"_biomarker_pairs.txt
			echo "done"
		else
			echo "file "$prefix"_sort.csv is empty. No Hits!"
			echo "Protein_ID,"$db","$db"_80/80" > "$3""$name"_best-hits.csv 
	fi		
			
	
	### remove file to save disk space and keep clarity
	rm "$prefix".csv
	rm "$prefix"_sort.csv
done

#### when all BLAST scores for all_biomarker are completed, merge them!
if [ -e "$3""$Venn_group"_vs_O26_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_O103_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O103-O121_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O104_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_O111_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_O111-O26_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O111-O26-O145_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O121_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O145_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_O157_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ST_11_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ST_11-ST_20-ST_29_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ST_20_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ST_20-ST_29_best-hits.csv ]  && [ -e "$3""$Venn_group"_vs_ST_29_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_CC_32_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_CC_32-ST_306_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_CC_32-ST_306-ST_678_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_CC_32-ST_678_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_Eae_pos_stx_neg_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_Genostar_biomarker_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_ST_306_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_ST_306-ST_678_best-hits.csv ] && [ -e "$3""$Venn_group"_vs_ST_678_best-hits.csv ]
	then 
		echo "------------all BLAST best hits done------------"
		echo "------------merge BLAST best hits---------------"
		### use a small Rscript to merge the five .csv files for each database together!
		Rscript "$4" "$3""$Venn_group"_vs_O26_best-hits.csv "$3""$Venn_group"_vs_O103_best-hits.csv "$3""$Venn_group"_vs_O103-O121_best-hits.csv "$3""$Venn_group"_vs_O104_best-hits.csv "$3""$Venn_group"_vs_O111_best-hits.csv "$3""$Venn_group"_vs_O111-O26_best-hits.csv "$3""$Venn_group"_vs_O111-O26-O145_best-hits.csv "$3""$Venn_group"_vs_O121_best-hits.csv "$3""$Venn_group"_vs_O145_best-hits.csv "$3""$Venn_group"_vs_O157_best-hits.csv "$3""$Venn_group"_vs_ST_11_best-hits.csv "$3""$Venn_group"_vs_ST_11-ST_20-ST_29_best-hits.csv "$3""$Venn_group"_vs_ST_20_best-hits.csv "$3""$Venn_group"_vs_ST_20-ST_29_best-hits.csv "$3""$Venn_group"_vs_ST_29_best-hits.csv "$3""$Venn_group"_vs_CC_32_best-hits.csv "$3""$Venn_group"_vs_CC_32-ST_306_best-hits.csv "$3""$Venn_group"_vs_CC_32-ST_306-ST_678_best-hits.csv "$3""$Venn_group"_vs_CC_32-ST_678_best-hits.csv "$3""$Venn_group"_vs_Eae_pos_stx_neg_best-hits.csv "$3""$Venn_group"_vs_Genostar_biomarker_best-hits.csv "$3""$Venn_group"_vs_ST_306_best-hits.csv "$3""$Venn_group"_vs_ST_306-ST_678_best-hits.csv "$3""$Venn_group"_vs_ST_678_best-hits.csv "$3""$Venn_group"_merged_vs_all.csv
		echo "------------"$Venn_group" BLAST results are merged-----------"
	else  
		echo "------------at least one BLAST best hit in "$Venn_group" still missing------------"
fi

	
