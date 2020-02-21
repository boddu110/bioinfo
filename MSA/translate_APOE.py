#!usr/bin/env python3
#translate_APOE.py

import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
protein=[]
for record in SeqIO.parse("APOE_refseq_transcript.fasta","fasta"):
    protein_seq=SeqRecord(record.seq.transcribe().translate(),id=record.id,description=record.description,name=record.name)
    protein.append(protein_seq)

SeqIO.write(protein,"apoe_aa.fasta","fasta")
