#!/usr/bin/env python3

import fileinput

# Write a program that reads a 23andme file
# Report the number of SNPs successfully assayed
# Report the heterozygosity
# Report the failure rate
# https://opensnp.org/genotypes (make sure the data is 23andme)

count = 0
fail = 0 # failure rate counter
hom = 0  # homozygosity counter
het = 0  # heterozygosity counter

for line in fileinput.input():
    if line.startswith ('#'): continue
    column = line.split()
    genotype = column[3]
    nt1 = genotype[0:1]
    nt2 = genotype[1:2]
    count += 1
    if   nt1 == '-':     fail += 1
    elif nt1 == nt2:     hom  += 1
    else:                het  += 1
snps = hom + het 
print(count, snps)   
print(f'%.2f %.2f' % (het/snps, fail/snps))


"""
proc23.py whatever.23andme.whatever
"""
