#!/usr/bin/env python3
#sliding_window.py
import sys
import re
from Bio import SeqIO
def sliding_window(k,string):
    kmers = []
    for start in range(0, len(string) - k + 1):
        kmers.append(string[start:start + k])
    return kmers

def gc_content(dna):
    '''returns [0,1] the fraction of gc's in given string '''
    dna = dna.lower()
    gc = 0
    for nucleotide in dna:
        if nucleotide in ['g','c']:
           gc += 1
        
    return gc/ len(dna)

if __name__ == "__main__":
      ''' find the length of sys arguments '''
      arg_count = len(sys.argv) -1
      if arg_count < 2:
          raise Exception("This script requires atleast 2 arguments")
      '''convert kmer size argument to int from a string '''
      k = int(sys.argv[1])
      fasta = sys.argv[2]

      for dna in SeqIO.parse(fasta, "fasta"):
          out = sliding_window(k, dna.seq)
      
      for kmer in out:
        print("{}\t{:.2f}".format(kmer,gc_content(kmer)))
