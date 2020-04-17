#!/usr/bin/env python3

# Write a program that prints the reverse-complement of a DNA sequence
# You must use a loop and conditional

dna = 'ACTGAAAAAAAAAAA'
comp = '' #complement
rcomp = '' #reverse complement

for nt in range(len(dna)):
    if   dna[nt] == 'A': comp += 'T'
    elif dna[nt] == 'T': comp += 'A'
    elif dna[nt] == 'C': comp += 'G'
    else:                comp += 'C'
print(dna)
print(comp)
print('---------------')

for nt in range(len(dna)):
    if   dna[nt] == 'A': rcomp = 'T' + rcomp
    elif dna[nt] == 'T': rcomp = 'A' + rcomp
    elif dna[nt] == 'C': rcomp = 'G' + rcomp
    else:                rcomp = 'C' + rcomp
print(dna)
print(rcomp)


"""
TTTTTTTTTTTCAGT
"""
