### Skript um die verschiedenen best Blast hits zu mergen!

library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
### load cammand line arguments as variables:

args <- commandArgs(trailingOnly = TRUE)
#print(args)
input1 <- args[1]
input2 <- args[2]
input3 <- args[3]
input4 <- args[4]
input5 <- args[5]
input6 <- args[6]
input7 <- args[7]
input8 <- args[8]
input9 <- args[9]
input10 <- args[10]
input11 <- args[11]
input12 <- args[12]
input13 <- args[13]
input14 <- args[14]
input15 <- args[15]
input16 <- args[16]
input17 <- args[17]
input18 <- args[18]
input19 <- args[19]
input20 <- args[20]
input21 <- args[21]
input22 <- args[22]
input23 <- args[23]
input24 <- args[24]

output_file <- args[25]

### read in input.csv files
O26 <- read.csv(input1, header = TRUE, sep = ",")
O103 <- read.csv(input2, header = TRUE, sep = ",")
O103_O121 <- read.csv(input3, header = TRUE, sep = ",")
O104 <- read.csv(input4, header = TRUE, sep = ",")
O111 <- read.csv(input5, header = TRUE, sep = ",")
O111_O26 <- read.csv(input6, header = TRUE, sep = ",")
O111_O26_O145 <- read.csv(input7, header = TRUE, sep = ",")
O121 <- read.csv(input8, header = TRUE, sep = ",")
O145 <- read.csv(input9, header = TRUE, sep = ",")
O157 <- read.csv(input10, header = TRUE, sep = ",")
ST11 <- read.csv(input11, header = TRUE, sep = ",")
ST11_ST20_ST29 <- read.csv(input12, header = TRUE, sep = ",")
ST20 <- read.csv(input13, header = TRUE, sep = ",")
ST20_ST29 <- read.csv(input14, header = TRUE, sep = ",")
ST29 <- read.csv(input15, header = TRUE, sep = ",")

CC32 <- read.csv(input16, header = TRUE, sep = ",")
CC32_ST306 <- read.csv(input17, header = TRUE, sep = ",")
CC32_ST306_ST678 <- read.csv(input18, header = TRUE, sep = ",")
CC32_ST678 <- read.csv(input19, header = TRUE, sep = ",")
IPEC_eae <- read.csv(input20, header = TRUE, sep = ",")
Genostar <- read.csv(input21, header = TRUE, sep = ",")
ST306 <- read.csv(input22, header = TRUE, sep = ",")
ST306_ST678 <- read.csv(input23, header = TRUE, sep = ",")
ST678 <- read.csv(input24, header = TRUE, sep = ",")
### merge files
output <- merge(ST11, O157, by="Protein_ID", all=TRUE)
output2 <- merge(output, ST20, by="Protein_ID", all=TRUE)
output3 <- merge(output2, ST29, by="Protein_ID", all=TRUE)
output4 <- merge(output3, ST20_ST29, by="Protein_ID", all=TRUE)
output5 <- merge(output4, ST11_ST20_ST29, by="Protein_ID", all=TRUE)
output6 <- merge(output5, O26, by="Protein_ID", all=TRUE)
output7 <- merge(output6, O103, by="Protein_ID", all=TRUE)
output8 <- merge(output7, O103_O121, by="Protein_ID", all=TRUE)
output9 <- merge(output8, O104, by="Protein_ID", all=TRUE)
output10 <- merge(output9, O111, by="Protein_ID", all=TRUE)
output11 <- merge(output10, O111_O26, by="Protein_ID", all=TRUE)
output12 <- merge(output11, O111_O26_O145, by="Protein_ID", all=TRUE)
output13 <- merge(output12, O121, by="Protein_ID", all=TRUE)
output14 <- merge(output13, O145, by="Protein_ID", all=TRUE)

output15 <- merge(output14, CC32, by="Protein_ID", all=TRUE)
output16 <- merge(output15, ST306, by="Protein_ID", all=TRUE)
output17 <- merge(output16, ST678, by="Protein_ID", all=TRUE)
output18 <- merge(output17, CC32_ST306, by="Protein_ID", all=TRUE)
output19 <- merge(output18, CC32_ST678, by="Protein_ID", all=TRUE)
output20 <- merge(output19, CC32_ST306_ST678, by="Protein_ID", all=TRUE)
output21 <- merge(output20, ST306_ST678, by="Protein_ID", all=TRUE)

output22 <- merge(output21, IPEC_eae, by="Protein_ID", all=TRUE)
output23 <- merge(output22, Genostar, by="Protein_ID", all=TRUE)

output23[is.na(output23)] <- 0

### rename header into a more human readable way:
names(output23) <- c("Protein_ID", "ST11", "ST11 80/80", "O157", "O157 80/80", "ST20", "ST20 80/80", "ST29", "ST29 80/80", "ST20_ST29", "ST20_ST29 80/80", "ST11_ST20_ST29", "ST11_ST20_ST29 80/80", "O26", "O26 80/80", "O103", "O103 80/80", "O103_O121", "O103_O121 80/80", "O104", "O104 80/80", "O111", "O111 80/80", "O111_O26", "O111_O26 80/80", "O111_O26_O145", "O111_O26_O145 80/80", "O121", "O121 80/80", "O145", "O145 80/80", "CC32", "CC32 80/80", "ST306", "ST306 80/80", "ST678", "ST678 80/80", "CC32_ST306", "CC32_ST306 80/80", "CC32_ST678", "CC32_ST678 80/80", "CC32_ST306_ST678", "CC32_ST306_ST678 80/80", "ST306_ST678", "ST306_ST678 80/80", "IPEC_eae", "IPEC_eae 80/80", "Genostar", "Genostar 80/80")

### write output file:
write.csv(output23, file = output_file, row.names = FALSE)
