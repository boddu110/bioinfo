#!/usr/bin/env python3
#sliding_window.py
import sys
   
def sliding_window(k,string):
    kmers=[]
    for start in range(0,len(string)-k+1):
        kmers.append(string[start:start+k])
    return kmers

def gc_content(dna):
    ''' Returns [0, 1], the fraction of GCs in the given string'''
    # For consistency, make the sequence lowercase
    dna = dna.lower()

    # Count the number of g's and c's
    gc = 0
    for nucleotide in dna:
        if nucleotide in ['g', 'c']:
            gc += 1

    return gc/len(dna)

if __name__ == "__main__":

    # Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires 2 arguments: a kmer size and then a string")

    for kmer in sliding_window(int(sys.argv[1]), sys.argv[2]):

        gc = round(gc_content(kmer),2)
        print("{}\t{}".format(kmer,gc))
