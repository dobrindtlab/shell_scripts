### Script to merge the given .csv files!

library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
### load cammand line arguments as variables:

args <- commandArgs(trailingOnly = TRUE)
#print(args)
input1 <- args[1]
input2 <- args[2]
input3 <- args[3]
input4 <- args[4]
input5 <- args[5]
output_file <- args[6]

### read in input.csv files
non_patho <- read.csv(input1, header = TRUE, sep = ",")
IPEC <- read.csv(input2, header = TRUE, sep = ",")
ExPEC <- read.csv(input3, header = TRUE, sep = ",")
IPEC_eae <- read.csv(input4, header = TRUE, sep = ",")
IPEC_ST11 <- read.csv(input5, header = TRUE, sep = ",")

### merge files
output <- merge(non_patho, ExPEC, by="Protein_ID", all=TRUE)
output2 <- merge(output, IPEC, by="Protein_ID", all=TRUE)
output3 <- merge(output2, IPEC_eae, by="Protein_ID", all=TRUE)
output4 <- merge(output3, IPEC_ST11, by="Protein_ID", all=TRUE)
output4[is.na(output4)] <- 0

### rename header into a more human readable way:
names(output4) <- c("Protein_ID", "non-patho", "non-patho 80/80", "ExPEC", "ExPEC 80/80", "IPEC", "IPEC 80/80", "IPEC_eae", "IPEC_eae 80/80", "IPEC_ST11", "IPEC_ST11 80/80")

### write output file:
write.csv(output4, file = output_file, row.names = FALSE)
