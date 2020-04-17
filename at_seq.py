#!/usr/bin/env python3

import random
#random.seed(1) # comment-out this line to change sequence each time

# Write a program that stores random DNA sequence in a string
# The sequence should be 30 nt long
# On average, the sequence should be 60% AT
# Calculate the actual AT fraction while generating the sequence
# Report the length, AT fraction, and sequence

length = 30
dna = ''
at_count = 0
for i in range(length):
    r = random.random()
    if r <= 0.6:
        r = random.randint(0,1)
        if r == 0: dna += 'A'
        else:      dna += 'T'
        at_count += 1
    else:
        r = random.randint(0,1)
        if r == 0: dna += 'C'
        else:      dna += 'G'
print(length, at_count / length, dna)

#alternative method of generating a random string of letters, but not AT rich
#alph = 'ATGC'
#for i in range(length):
#    nt = random.choice(alph)
#    dna += nt
#    at_count += 1
#print(dna)



"""
30 0.6666666666666666 ATTACCGTAATCTACTATTAAGTCACAACC
"""
