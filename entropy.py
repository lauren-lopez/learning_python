#!/usr/bin/env python3

# Write a Shannon entropy calculator: H = -sum(pi * log(pi))
# Use fileinput to get the data from nucleotides.txt
# Make sure that the values are probabilities
# Make sure that the distribution sums to 1
# Report with 3 decimal figures

import fileinput
import math

# Need to split the data --> currently looks like 'A 0.1' etc on each line. We only want the '0.1' part.
# Need to process the part to grab the part that we want out of it.

data = [] 
sum = 0                              # creates an array to put the data into
for line in fileinput.input():          # for each line in the file we are working with via the command line
    if line.startswith('#'): continue   # ignores any comments starting w/ #
    values = line.split()               # splits the strings into columns (splits on spaces) --> assigns each word/number etc to 'values'
    v = float(values[1])                # v = whatever 'values' contains at position 1 (remember, starts at 0) and makes it a float
    assert(v >= 0 and v <= 1)           # ensure 0>v<1
    sum += v                            # sum all 'v'
    data.append(v)                      # appends all 'v' into a list
assert(math.isclose(sum, 1))            # ensures the sum is close to 1
print(data)

# Shannon Entropy

h = 0
for i in range(len(data)):
    #h += data[i] * math.log2(data[i])   # entropy = sum of (data at position[i] * log2 data at [i])
    h -= data[i] * math.log2(data[i])    # output is negative, so to make positive, either subtract all values (like here)
print(f'%.3f' % (h))                                 # or just multiply by -1



"""
python3 entropy.py nucleotides.txt
1.846
"""
