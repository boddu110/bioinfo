-!/usr/bin/env python3
wBioPython_seqio.py
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
if __name__ == "__main__":
	input=sys.argv[1]
	output=sys.argv[2]
#creating the list
revseq=[]
#reading the fasta file using SeqIO object
for rec in SeqIO.parse(input, "fasta"):
	sequ=rec.seq
        #using the reverse complement function
	sequ=sequ.reverse_complement()
        #updating sequ 
	rec.seq=sequ
        #updating all the values to the list
	revseq.append(rec)
SeqIO.write(revseq, output, "fasta")

