import sys
from Bio import SeqIO

fasta_input = open(sys.argv[1], "r")
fasta_output = open(sys.argv[2], "a")

fasta_ID = sys.argv[3]
fasta_start = int(sys.argv[4])
fasta_stop = int(sys.argv[5])


with fasta_output as f:
	for seq_record in SeqIO.parse(fasta_input, "fasta"):
		if seq_record.id == fasta_ID and int(fasta_start)<int(fasta_stop):
			f.write(">"+str(seq_record.id)+"_uidA"+"\n")
			fasta_start_1 = fasta_start - 1
			f.write(str(seq_record.seq[int(fasta_start_1):int(fasta_stop)]) + "\n")  #write specified subsequence
			
		elif seq_record.id == fasta_ID and int(fasta_start)>int(fasta_stop):
			f.write(">"+str(seq_record.id)+"_uidA"+"\n")
			fasta_stop_1 = fasta_stop - 1
			DNA=seq_record.seq[int(fasta_stop_1):int(fasta_start)]
			f.write(str(DNA.reverse_complement()) + "\n")  #if the orientation of uidA is false
			
