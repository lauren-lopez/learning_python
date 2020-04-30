#!/usr/bin/env python3

import fileinput

# Write a program that computes typical sequence stats
# No, you cannot import any other modules!
# Use rand_seq to generate the sequences
# Expected output is shown below

line_count, sum_nt, A_count, T_count, G_count, C_count = 0, 0, 0, 0, 0, 0

nt_lengths = []
for line in fileinput.input():
    if line.startswith('>'): continue   # removes '>seq' headers for each sequence  
    line = line.rstrip()                # removes return characters  
    line_count += 1
    nt = len(line)
    sum_nt += nt
    nt_lengths.append(nt)

    nt_list = []                             # splits sequence in a string into individual
    nt_list[:0] = line                       # nucleotides in a list

    for c in nt_list:                        # counts up numbers of each nucleotide in all sequences
        if   c == 'A': A_count += 1
        elif c == 'T': T_count += 1
        elif c == 'G': G_count += 1
        else:          C_count += 1

A_frac = A_count/sum_nt
T_frac = T_count/sum_nt
G_frac = G_count/sum_nt
C_frac = C_count/sum_nt

nt_lengths.sort(reverse = True)             # sorts lengths of sequences in decending order
min = nt_lengths[-1]
max = nt_lengths[0]

n50_sum = 0
for l in nt_lengths:                        # for each length in list of lengths,
    n50 = False
    n50_sum += l                            # summing lengths from largest to smallest
    if n50_sum > (sum_nt/2):                # running sum must equal half the total sum                           
        n50 = True
        break
#    if n50:
#        n50_len = l     # need to select the shortest length after the sum of lengths equals half of the total lengths


# Output
print(f'Number of sequences: {line_count}')
print(f'Number of letters: {sum_nt}')
print(f'Minimum length: {min}')
print(f'Maximum length: {max}')
#print(f'N50: {n50_len}')
print(f'Composition: A={A_frac:.3f} C={C_frac:.3f} G={G_frac:.3f} T={T_frac:.3f}')


"""
python3 rand_seq.py 100 100 100000 0.35 | python3 seqstats.py
Number of sequences: 100
Number of letters: 4957689
Minimum length: 219
Maximum length: 99853
N50: 67081
Composition: A=0.325 C=0.175 G=0.175 T=0.325
"""
