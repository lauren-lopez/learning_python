#!/usr/bin/env python3

import gzip
import sys
import biotools as bt
import argparse

# Use argparse
# Compute gc and gc-skew (g-c)/(g+c)

# setup
parser = argparse.ArgumentParser(
	description='Calculates GC skew within a window of a DNA sequence.')
# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='input DNA file')
parser.add_argument('--win', required=True, type=int,
	metavar='<int>', help='window size')
# finalization
arg = parser.parse_args()

for name, seq in bt.read_fasta(arg.file):
    for i in range(len(seq) -arg.win+1):
        sseq = seq[i:i+arg.win]
        print(f'{name}\t{i}\t{bt.gc(sseq):.2f}\t{bt.gc_skew(sseq):.3f}')

"""
python3 skew.py --file genome.fa.gz --win 100 | head
I	0	0.510	-0.333
I	1	0.500	-0.360
I	2	0.490	-0.347
I	3	0.490	-0.306
I	4	0.500	-0.320
I	5	0.510	-0.333
I	6	0.510	-0.333
I	7	0.500	-0.360
I	8	0.490	-0.347
I	9	0.490	-0.306
"""
