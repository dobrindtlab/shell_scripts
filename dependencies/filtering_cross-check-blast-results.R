### Skript zur bearbeitung von BLAST ergebnissen, um den besten Hit zu bekommen. 
library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")

### einlesen der Datei:
### Dabei darauf achten, schon vorher den BLAST Output in Excel leicht zu bearbeiten.
### BLAST output sortieren, bei identity "." durch "," ersetzen
### und bei queryacc den :(number) Anhang entfernen.
### als .csv speichern, but ";" als field seperator --> can man ändern!

args <- commandArgs(trailingOnly = TRUE)
#print(args)
input <- args[1]
venn_name <- args[2]
venn_name_80 <- args[3]
db <- args[4]
db_80 <- args[5]
input_strains <- args[6]
output_file <- args[7]
biomarker <- args[8]

blast <- read.csv(input, header = FALSE, sep = ";")

### rename the columns in a more readable format:

names(blast) <- c("queryid", "query", "strainid", "identity", "alignment_length", 
                  "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end",
                  "evalue", "bitscore", "coverage") 

### extract the strain names from the strainid
strains <- sub("(_[0-9]+)$", "", blast$strainid)

### add the column strains to the data.frame
blast2 <- mutate(blast, strain = strains)

### concentrate on the really important colums:

sumblast <- select(blast2, queryid, query, strainid, identity, bitscore, coverage, strain)

### group the data_frame on basis of equel strain IDs:

strain_group <- group_by(sumblast, query, strainid)

### display only the rows with the highest bitscore = highest identity
### additionally removes all hits with ident/cov lower than 80%
blast_hits_all <- filter(strain_group, bitscore == max(bitscore))
blast_hits_80_80 <- filter(strain_group, bitscore == max(bitscore), coverage > 80, identity > 80)

### load in protein_IDs for the given Venn group:

strain_IDs <- read.table(input_strains, header = TRUE, sep = ";")
strain_name <- as.data.frame(strain_IDs$strain)
names(strain_name) <- c("strain")

### merge all Blast hits to the protein_IDs
venn_name <- merge(strain_IDs, data.frame(table(strain = blast_hits_all$strain)), all = TRUE)
venn_name_80 <- merge(strain_name, data.frame(table(strain = blast_hits_80_80$strain)), all = TRUE)

### merge the two dataframes:
output <- merge(venn_name, venn_name_80, by=1, all=TRUE)

### rename the column names
names(output) <- c("hit in strain", "Venn-Group", "Venn-Group 2", db, db_80)
output[is.na(output)] <- 0
### done!

write.csv(output, file = output_file, row.names = FALSE)

### get biomarker pairs, to identify duplicate biomarker
biomarker_pair <- select(blast_hits_80_80, query, strain, bitscore)

write.table(biomarker_pair, file = biomarker, row.names = FALSE, sep = "\t")
