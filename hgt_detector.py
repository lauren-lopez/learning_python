#!/usr/bin/env python3

import argparse
import biotools
import re
import math
import json

parser = argparse.ArgumentParser(
    description='Horiztonal gene transfer detector.')
parser.add_argument('--file', required=True, type=str,
    metavar='<str>', help='FASTA file')
parser.add_argument('--pseudo', required=False, type=float, default=1.0,
    metavar='<float>', help='pseudocount [%(default)f]')
arg = parser.parse_args()


dd = {}                                                 # Dictionary is 2D and contains codon counts for each locus
genome = {}                                             # Dictionary is 1D and contains total codon counts for whole genome
sizes = {}                                              # Number of codons in sequence at each locus    
total = 0                                               # Total counts of all codons in the genome
for name, seq in biotools.read_fasta(arg.file):
    match = re.search('locus_tag=(\w+)', name)          # Find locus tag
    locus = match[1]
    sizes[locus] = 0
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if locus not in dd:                             # Adds locus name to dictionary
            dd[locus] = {}
        if codon not in dd[locus]:                      # First time a codon is seen at a locus, need to add pseudocount to both codon count and size count
            dd[locus][codon] = arg.pseudo + 1               # Adds pseudocount plus first codon count
            sizes[locus] += arg.pseudo + 1                  # Increases size length by pseudocount plus first codon count
        else:                                           # Every other time codon is seen at a locus, increase codon count and size count by one
            dd[locus][codon] += 1
            sizes[locus] += 1
        
        if codon not in genome:                         # First time a codon is seen, need to add baseline pseudocount for each codon, then from there sum up 
            genome[codon] = arg.pseudo + 1                  # Counts for each codon
            total += arg.pseudo + 1                         # Total codons in genome
        else:                   
            genome[codon] += 1
            total += 1


sum = 0   
for codon in genome:                                # For each codon in total genome dictionary
    for locus in dd.keys():                         # For each locus in the dictionary that stores codon counts at that locus
        if codon not in dd[locus]:                  # If the genome codon is not seen at that locus, 
            dd[locus][codon] = arg.pseudo               # then need to store it as a pseudocount
            sizes[locus] += arg.pseudo                  # and also need to add pseudocount to codon length
    genome[codon] /= total                          # Converts codon counts within genome to probabilities
    sum += genome[codon]
assert(math.isclose(sum, 1))                        # All probabilities should sum up to 1


for locus in dd.keys():                             # For each locus in the dictionary that stores codon counts at that locus
    sum = 0
    for codon in dd[locus]:                         # For each codon at locus
        dd[locus][codon] /= sizes[locus]            # Converts codon counts at locus to probabilities
        sum += dd[locus][codon]
    assert(math.isclose(sum, 1))                    # All probabilities should sum up to 1


kl_distance = {}
for locus in dd.keys():                             # For each locus in the dictionary that stores codon probabilities at that locus
    kl = 0   
    for codon in genome.keys():                     # For each codon in the genome dictionary
        kl += genome[codon] * math.log2((genome[codon]/dd[locus][codon]))   # Comparing genome to gene
    kl_distance[locus] = kl                         # Dictionary stores the kl distance at that locus    

kl_distance = sorted(kl_distance.items(), key = lambda x: x[1], reverse = True)

for kld in kl_distance:    
    print(f'Locus ID:{kld[0]}\tKL Distance:{kld[1]:.3f}')

    



# for each locus, take the codon count and divide by number of possible codons (len(seq)/3)










    #pseudocount - add 1 to all codon counts (like a baseline)
    #codon probabilities - divide by num codons in each gene individually
    #codon probabilities for entire genome - divide by all codons seen
    #P = codon probability for whole genome
    #Q = codon probaility for each gene
