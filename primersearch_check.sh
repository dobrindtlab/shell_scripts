#!/bin/bash
set -e

### script to test primer for biomarker against different test genomes
### use primersearch and check the output automatically

if [ $# -lt 3 ]; then
	echo "Usage: $0 input_genomes primer_folder output_folder"
	exit 1
fi

### primersearch primer_pair infile needs to be a 3 column tab-seperated txt file with name, forward and reverse primer.
for file in "$1"*.fna; do
	prefix="${file%*.fna}"
	name="${prefix##*/}"
	
	for primer in "$2"*_primer.txt; do
		pre="${primer%*_primer.txt}"
		primer_name="${pre##*/}"
		echo "for genome "$name" and primer "$primer_name""
		primersearch -seqall "$file" -infile "$primer" -mismatchpercent 0 -outfile "$3""$name"_vs_"$primer_name"_primer.out
	
	done
done

### generate an output which is easier to read on a quick glance!
### run primersearch, grep only lines 1 above and 5 below "Amplimer 1" 

echo "clean primersearch result"
for file in "$3"*primer.out; do
	prefix="${file%*.out}"
	name="${prefix##*/}"
	primer="${prefix##*_vs_}"
	[ -d "$3""$primer"/ ] ||  mkdir "$3""$primer"/
	if grep -q 'Amplimer 1' < "$file"; then
		
		grep -A 5 -B 1 'Amplimer 1' < "$file" > "$3""$primer"/"$name"_clean.out
		
	else
		echo "no Primer hits found" > "$3""$primer"/"$name"_clean.out
	fi
	rm "$file"
	 
done
