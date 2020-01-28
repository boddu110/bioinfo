#!/usr/bin/env python3
#BioPython_seq.py
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


seq_s=Seq("aaaatgggggggggggccccgtt",generic_dna)
sequence_r=SeqRecord(seq_s,id="12345",description="example1")
print(sequence_r)
#using SeqIO object to write a file
SeqIO.write(sequence_r,"BioPython_seq.gb", "gb")		
