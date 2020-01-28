#!/usr/bin/env python
#BioPython_genbank.py
from Bio import SeqIO
from Bio import Entrez

ntrez.email = "boddu.j@northeastern.edu"
seq=[]
#fetching the sequence using Entrez.efetch obtaining the genbank entry with id 515056
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="515056") as handle:
	seq_record = SeqIO.read(handle, "gb")
	seq.append(seq_record)
#fetching the sequence using Entrez.efetch using accession id J01673.1
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="J01673.1") as handle:
	seq_record = SeqIO.read(handle, "gb") 
	seq.append(seq_record)
print(seq)

for i in seq:
	for entry in i.features:
		print("{}\t{}\t{}".format(entry.type,entry.location,entry.strand))

