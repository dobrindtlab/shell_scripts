#!/bin/bash
set -e

### Identify biomarker after Prokka run:
### combine all different scripts into one.
### all depended scripts need to be in the same folder

if [ $# -lt 1 ]; then
	echo "Usage: $0 input-folder"
	exit 1
fi

project_name="${1%*/}"
name="${project_name##*/}"
	### löscht alle Zeichen bis zum letzen /

### STEP 1:
### extract protein seq from Venn groups: /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/extract_proteins_seq_of-Venn_groups.sh
### Usage: $0 input-folder project-name output_folder[absolut path]

bash extract_proteins_seq_of-Venn_groups.sh "$1" "$name" "$1"extracted-proteins_strict_80-40


### STEP 2:
### run BLAST on these proteins: /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/blast_Venn-groups_vs_control-DBs.sh
### Usage: $0 input_folder[extracted-proteins] inputfolder[BLAST DB] output-path[Blast]

bash blast_Venn-groups_vs_control-DBs.sh "$1"extracted-proteins_strict_80-40 BLAST "$1"BLAST/

### STEP 3:
### run Rscript to analyse BLAST results: /mnt/WD2/matthias/own_Genostar_NCBI_NVI/Prokka/filtering_blast-results.sh
### Usage: $0 path-to-filtering-R-script path-to-BLAST-results output-path[Blast] path-to-merge-R-script path-to-Protein-IDs

bash filtering_blast-results.sh filtering_blast-results.R "$1"BLAST/ "$1"BLAST/results/ merge_blast-results.R "$1"extracted-proteins_strict_80-40/

echo "This Job is done. AWESOME!!!"


