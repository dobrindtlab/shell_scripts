#!/bin/bash
set -e

### script to gather specific proteins from Venn diagram groups and paste them in a new file

if [ $# -lt 2 ]; then
	echo "Usage: $0 input-folder project-name output_folder[absolut path]"
	exit 1
fi

[ -d "$3" ] ||  mkdir "$3"

if [ -e "$3"/"all_proteins.faa" ]  
	then 
		echo "------------protein index already exists------------"
	else  
		echo "create protein index for all Prokka hits"
		cat "$1"output_Prokka/*.faa > "$3"/all_proteins.fasta
		echo "done"
fi


### I want to analyse all Venn groups which are not shared with the non-patho control group.
### The Perl Script po2group_stats.pl names all .tsv files according to the group file
### extract the group_file header to use them to name the correct files.

group_file="$1"group_file.tsv
c=$(head -1 "$group_file" | tr '\t' '\n' | wc -l)

echo "------------ collect group_file header ------------"
if (( "$c" == 2 ))
	then
		g1=$(cut -f1 "$1"group_file.tsv | head -n 1)
		echo "------------extract protein IDs from specific files------------"
		cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"_protein-IDs.txt
		echo "------------ done ------------"
		echo "------------extract proteins from protein index with samtools------------"
		xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"_protein-IDs.txt > "$3"/"$g1"_proteins.faa
		echo "------------ done ------------"
	else 
		if (("$c" == 3 ))
			then
				g1=$(cut -f1 "$1"group_file.tsv | head -n 1)
				g2=$(cut -f2 "$1"group_file.tsv | head -n 1)
				echo "------------extract protein IDs from specific files------------"
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"_protein-IDs.txt
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g2"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g2"_protein-IDs.txt 
				echo "------------ done ------------"
				echo "------------extract protein IDs from two-fold shared files------------"
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"-"$g2"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"-"$g2"_protein-IDs.txt 
				echo "------------ done ------------"
				echo "------------extract proteins from protein index with samtools------------"
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"_protein-IDs.txt > "$3"/"$g1"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g2"_protein-IDs.txt > "$3"/"$g2"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"-"$g2"_protein-IDs.txt > "$3"/"$g1"-"$g2"_proteins.faa
			else
					
				g1=$(cut -f1 "$1"group_file.tsv | head -n 1)
				g2=$(cut -f2 "$1"group_file.tsv | head -n 1)
				g3=$(cut -f3 "$1"group_file.tsv | head -n 1)
				g4=$(cut -f4 "$1"group_file.tsv | head -n 1)
				echo "------------ done ------------"

### cut takes the second tab seperated column which is the protein ID
### tail -n +2 is needed to delete the header line
#### I need to extract the g1-3 specific files shared g1/2 g2/3 g3/1 not-in-g4

				echo "------------extract protein IDs from specific files------------"

				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"_protein-IDs.txt 
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g2"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g2"_protein-IDs.txt 
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g3"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g3"_protein-IDs.txt 

				echo "------------ done ------------"

				echo "------------extract protein IDs from two-fold shared files------------"

				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"-"$g2"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"-"$g2"_protein-IDs.txt 
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g1"-"$g3"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g1"-"$g3"_protein-IDs.txt 
				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g2"-"$g3"_specific_OGs.tsv |  tail -n +2 > "$3"/"$g2"-"$g3"_protein-IDs.txt 

				echo "------------ done ------------"

				echo "------------extract protein IDs from three-fold shared files------------"

				cut -f2 "$1"ORF_results/"$2"_i80_c40/strict_100/"$g4"_absent_OGs.tsv |  tail -n +2 > "$3"/"$g1"-"$g2"-"$g3"_protein-IDs.txt 

				echo "------------ done ------------"

### with samtools extract aa sequences of proteins using protein-IDs.txt files.

				echo "------------extract proteins from protein index with samtools------------"
				
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"_protein-IDs.txt > "$3"/"$g1"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g2"_protein-IDs.txt > "$3"/"$g2"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g3"_protein-IDs.txt > "$3"/"$g3"_proteins.faa
				
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"-"$g2"_protein-IDs.txt > "$3"/"$g1"-"$g2"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"-"$g3"_protein-IDs.txt > "$3"/"$g1"-"$g3"_proteins.faa
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g2"-"$g3"_protein-IDs.txt > "$3"/"$g2"-"$g3"_proteins.faa
				
				xargs samtools faidx "$3"/all_proteins.fasta < "$3"/"$g1"-"$g2"-"$g3"_protein-IDs.txt > "$3"/"$g1"-"$g2"-"$g3"_proteins.faa
				echo "------------ done ------------"
		fi
fi		
