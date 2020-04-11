#!/usr/bin/env python3

# Print out all the codons for the sequence below in reading frame 1
# Use a 'for' loop

dna = 'ATAGCGAATATCTCTCATGAGAGGGAA'

for i in range(0, len(dna), 3):
    codon = dna[i:i+3]    
    print(codon)
 
""" 

#Additional task: print out the codons for reading frames 1, 2, and 3, with labels

for f in range(3):
    print('reading frame', f+1)
    for i in range(f, len(dna) -2,3):
        codon = dna[i:i+3]
        print(codon)
        
#Additional task: print out all the kmers for k = 10

k = 10
for i in range(0, len(dna) -k+1):
    kmer = dna[i:i+k]
    print(kmer)





ATA
GCG
AAT
ATC
TCT
CAT
GAG
AGG
GAA
"""
