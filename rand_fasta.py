#!/usr/bin/env python3

import gzip
import sys
import math
import random

# Write a program that creates random fasta files
# Create a function that makes random DNA sequences
# Parameters include length and frequencies for A, C, G, T
# Command line: python3 rand_fasta.py <count> <min> <max> <a> <c> <g> <t>

assert len(sys.argv) == 8

def rand_dna(lenseq, a, c, g, t):
    dnaseq = []
    for i in range(lenseq):
        r = random.random()
        if   r < a:     dnaseq.append('A')
        elif r < a+c:   dnaseq.append('C')
        elif r < a+c+g: dnaseq.append('G')
        else:           dnaseq.append('T')
    return(''.join(dnaseq))
    
count = int(sys.argv[1])
min = int(sys.argv[2])
max = int(sys.argv[3])
a = float(sys.argv[4])
c = float(sys.argv[5])
g = float(sys.argv[6])
t = float(sys.argv[7])

for i in range(count):
    l = random.randint(min, max)
    dna = rand_dna(l, a, c, g, t)
    print(f'>Seq-{i+1}')
    print(dna)

