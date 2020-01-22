#!/usr/bin/env python
#hamming.py
import sys
def hamming(string1,string2):
    '''Returns hamming distance of two sequences of equal length'''
    z = 0
    for a, b in zip(string1, string2):
        if a != b:
            z+= 1
    return z
if __name__ == "__main__":
     arg_count = len(sys.argv) -1
     if arg_count <2:
         raise Exception("This script requires two argument for equal length")
     x = sys.argv[1]
     y = sys.argv[2]
     if len(x) != len(y):
         raise ValueError("Sequence need to be the same length")
     result = hamming(x,y)
     print('{}\t{}\t{}'.format(x,y,result))
