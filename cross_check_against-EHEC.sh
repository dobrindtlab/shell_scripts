#!/bin/bash
set -e

###	check biomarker against all other strains!
###	Find out if biomarker are present in other strains, and if so in which Venn group they are!
### probably combination of samtools, cut, BLASTp and Rstudio scripts
### need to exclude own group, may be build a database for each venngroup and merge results
### need to know the number of genomes in each venngroup if present only in 1/20 it is still ok.

if [ $# -lt 2 ]; then
	echo "Usage: $0 input_path filtering_cross-check-blast-results.R_path"
	exit 1
fi

####### 	Notes for use		######
### Need to create the strains.txt files in advance! Also all shared Venn groups!
### create the biomarker.index file for every group!

for file in "$1"*strains.txt; do
	prefix="${file%*_strains.txt}"
	Venn_group="${prefix##*/}"
	
	[ -d "$1"all-stx_w-o_"$Venn_group"/ ] ||  mkdir "$1"all-stx_w-o_"$Venn_group"/
	### By default, comm outputs 3 columns: left-only, right-only, both. The -1, -2 and -3 switches suppress these columns.
	### e.g. comm -23 <(sort -u all-stx_strains.txt) <(sort -u ST11_strains.txt) > all-stx_w-o_ST11.txt
	### it is important to sort the text files before!
	comm -23 <(sort -u "$1"all-stx.txt) <(sort -u "$file") > "$1"all-stx_w-o_"$Venn_group"/all-stx_w-o_"$Venn_group".index

	for file in "$1"all-stx_w-o_"$Venn_group"/*.index; do
		### take the newly created index file to copy all proteins of each strain into the new folder
		echo "------ copy proteins for "$Venn_group" -------"
		while read -r line; do
		cp /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/all_250_strains/output_Prokka/"$line".faa "$1"all-stx_w-o_"$Venn_group"/"$line".faa;
		done < "$file"
	done
	cat "$1"all-stx_w-o_"$Venn_group"/*.faa > "$1"all-stx_w-o_"$Venn_group"/all-stx_w-o_"$Venn_group".fasta
	rm "$1"all-stx_w-o_"$Venn_group"/*.faa
	
	### build up the BLAST DB for each Venn_group:
	echo "------ make BLAST DB for "$Venn_group" -------"
	makeblastdb -in "$1"all-stx_w-o_"$Venn_group"/all-stx_w-o_"$Venn_group".fasta -dbtype 'prot' -title all-stx_w-o_"$Venn_group" -out "$1"all-stx_w-o_"$Venn_group"/all-stx_w-o_"$Venn_group"
done

### copy the biomarker which are of interest for each group and put them in the same folder
### need to extract the biomarker from former biomarker files!

for file in "$1"*.index; do
	prefix="${file%*_biomarker.index}"
	name="${prefix##*/}"
	
	echo "------ extracting biomarker for "$name" -------"
	xargs samtools faidx /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/cross_check_biomarker/"$name"_biomarker.faa < "$file" > "$1""$name"_biomarker.faa
	
done


### now BLASTp all biomarker proteins of the Venn_Group against the other proteomes of all other strains!

for dir in "$1"all-stx_w-o*/; do
	dir_name="${dir%*/}"
	dir_name="${dir_name##*/}"
	venn_group="${dir_name##*w-o_}"
	echo "------- BLASTp of "$venn_group" biomarker against DB "$dir_name""
	[ -d "$1"BLAST_results/ ] ||  mkdir "$1"BLAST_results/
	blastp -task blastp -query "$1""$venn_group"_biomarker.faa -db "$dir"/"$dir_name" -evalue 1e-9 -num_threads 32 -parse_deflines -outfmt '7 qseqid qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out "$1"BLAST_results/"$venn_group"_vs_"$dir_name".out

done


### now filter the results  
[ -d "$1"BLAST_results/filtered/ ] ||  mkdir "$1"BLAST_results/filtered/
for file in "$1"BLAST_results/*.out; do
	prefix="${file%*.out}"
	### behält alle Zeichen vor dem .out
	name="${prefix##*/}"
	### löscht alle Zeichen bis zum letzen /
	db="${name##*vs_}"
	Venn_group="${name%*_vs*}"
	
	### change from tab delimited file to ";" delimited file:
	sed 's/\t/;/g' "$file" > "$prefix".csv
	
	### delete all line with "#" at the start --> not important information and could not read by R
	sed '/^#/ d' "$prefix".csv > "$prefix"_sort.csv
	
	### input prepared file into R for analysis:
	### usage: R R_schript BLAST_results.csv Name_der_BLAST_kombi output-PATH
	### Problem: If no hits were found in the VENN analysis, files are empty! Problem with R!
	
	if [[ -s "$prefix"_sort.csv ]] ### checks if file is empty or not
		then
			echo "run Rscript for "$Venn_group" with database "$db"."	
			Rscript "$2" "$prefix"_sort.csv "$name" "$name"_80_80 "$db" "$db"_80/80 "$1"all-stx_strains.csv "$1"BLAST_results/filtered/"$name"_best-hits.csv "$1"BLAST_results/filtered/"$name"_query-strain_pair.txt
			echo "done"
		else
			echo "file "$prefix"_sort.csv is empty. No Hits!"
			echo "Protein_ID,"$db","$db"_80/80" > "$1"BLAST_results/filtered/"$name"_best-hits.csv 
	fi		
			
	
	### remove file to save disk space and keep clarity
	rm "$prefix".csv
	rm "$prefix"_sort.csv
done







