#!/usr/bin/env python3

# Write a program that computes the GC% of a DNA sequence
# Format the output for 2 decimal places
# Use all three formatting methods

dna = 'ACAGAGCCAGCAGATATACAGCAGATACTAT' # feel free to change

#Code written independently prior to class
GC_count = 0 #counting number of G and C in dna

for nt in range(len(dna)):
    if   dna[nt] == 'G': GC_count += 1
    elif dna[nt] == 'C': GC_count += 1             
    GC_content = GC_count/len(dna)

print('formatting method 1')
print('%.2f' % (GC_content))
print('formatting method 2')
print('{:.2f}' .format(GC_content))
print('formatting method 3')
print(f'{GC_content:.2f}')


#Code written in class
#gc_count = 0
#for c in dna:
#    if c == 'G' or c == 'C':
#        gc_count += 1
#print(gc_count/len(dna))
#print('%.2f' % (gc_count/len(dna)))
#print('{:.2f}' .format(gc_count/len(dna))
#print(f'{(gc_count/len(dna)):.2f}')

"""
0.42
0.42
0.42
"""