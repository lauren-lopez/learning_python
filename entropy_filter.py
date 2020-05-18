#!/usr/bin/env python3

# Write a program that masks areas of low complexity sequence
# Use argparse for command line arguments (see example below)
# Use read_fasta() from biotools.py

import argparse
import biotools
import math

# setup
parser = argparse.ArgumentParser(
	description='Brief description of program.')
# required arguments
parser.add_argument('--input', required=True, type=str,
	metavar='<path>', help='input DNA sequence (FASTA file)')

# optional arguments with default parameters
parser.add_argument('--win', required=False, type=int, default=15,
	metavar='<int>', help='window size [%(default)i]')
parser.add_argument('--threshold', required=False, type=float, default=1.1,
	metavar='<float>', help='entropy threshold [%(default)f]')

# switches
parser.add_argument('--lowercase', action='store_true',
	help='report lowercase instead of N')
# finalization
arg = parser.parse_args()

def entropy(data):
    h = 0
    for i in range(len(data)):                  # data is list of fractions of each nucleotide
        h -= data[i] * math.log2(data[i])        
    return(h)

for name, seq in biotools.read_fasta(arg.input):
    filtered_seq = list(seq)
    for i in range(0, len(seq) -arg.win+1):
        sseq = seq[i : i+arg.win]
        a, t, g, c = 0.0, 0.0, 0.0, 0.0
        for nt in sseq:
            if   nt == 'A': a += 1.0
            elif nt == 'T': t += 1.0
            elif nt == 'G': g += 1.0
            elif nt == 'C': c += 1.0
            
        a_frac = a/arg.win
        t_frac = t/arg.win
        g_frac = g/arg.win
        c_frac = c/arg.win
            
        data = []                               # creates a list of the fractions of nucleotides within a window
            
        if a_frac > 0: data.append(a_frac)      # only appends fractions that are > 0 (otherwise the entropy calculator does not work)
        if t_frac > 0: data.append(t_frac)
        if g_frac > 0: data.append(g_frac)
        if c_frac > 0: data.append(c_frac)
        
        sseq_entropy = float(entropy(data))     # calculates entropy of window
            
        
        if sseq_entropy < arg.threshold:                            # if entropy in window is less than threshold,
            if arg.lowercase:                                       # change each character within that window 
                for j in range(i, i+arg.win):                       #   to either lowercase letter or 'N'
                    filtered_seq[j] = filtered_seq[j].lower()       # Note we are changing the LIST version of the sequence
            else:
                for j in range(i, i+arg.win):
                    filtered_seq[j] = 'N'
                
    print(f'>{name}')
    print(''.join(filtered_seq))         
            
  




"""
python3 entropy_filter.py --help
usage: entropy_filter.py [-h] --input <path> [--window <int>]
                         [--threshold <float>] [--lowercase]

Low complexity sequence masker.

optional arguments:
  -h, --help           show this help message and exit
  --input <path>       fasta file
  --window <int>       optional integer argument [15]
  --threshold <float>  entropy threshold [1.100000]
  --lowercase          report lower case instead of N


python3 entropy_filter.py --input genome.fa.gz | head -20
>I
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAAAAATTGAGATAAGAAAACATTTTACTTTTTCAAAATTGTTTTCATGC
TAAATTCAAAACNNNNNNNNNNNNNNNAAGCTTCTAGATATTTGGCGGGTACCTCTAATT
TTGCCTGCCTGCCAACCTATATGCTCCTGTGTTTAGGCCTAATACTAAGCCTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAAAAGAATATGGTAGCTACAGAAACGGTAGTACACTCTTCTGNNNNNNNNNNNNNN
NTGCAATTTTTATAGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAAT
TGCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAAACANNNNNNNNNNNNNNNGAAAT
GAATATCGTAGCTACAGAAACGGTTGTGCACTCATCTGAAANNNNNNNNNNNNNNNNNNN
NNGCACTTTGTGCAGAATTCTTGATTCTTGATTCTTGCAGAAATTTGCAAGAAAATTCGC
"""
