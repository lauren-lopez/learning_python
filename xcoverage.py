#!/usr/bin/env python3

# Write a program that simulates random BAC coverage over a genome
# Command line arguments include
# 	Genome size (e.g. 1000)
# 	X coverage (e.g. 5)
# Use assert() to check parameter bounds
# Report min, max, and histogram of coverage
# Note that your output may vary due to random function

# randomly sampling the genome, and hoping to cover each position in the genome
# 5x coverage for a genome size of 1000 means 5000 reads
# some positions will be samples multiple times, some positions will not be sampled at all

import sys
import random

assert(len(sys.argv) == 3)

gsize = int(sys.argv[1])
x = float(sys.argv[2])

assert(gsize > 0)
assert(x > 0)

bacs = int(gsize * x) # the coverage of the genome (number of samples)
genome = [0] * gsize  # this creates an array of positions in the genome that can be sampled
for i in range(bacs): # this sets the number of samples to be taken (genome coverage)
    r = random.randint(0, gsize-1) # generates a random number within the size of the genome
    genome[r] += 1    # the random poisition within the genome (generated above) gains +1 each time it is sampled

genome.sort()    # sort the numbers of reads at each position in the genome
min = genome[0]  # position with the fewest number of reads
max = genome[-1] # position with the most number of reads will be at the end [-1] of the sort 

hist = [0] * (max + 1) # generates an array for the histogram
for v in genome:       # for value in genome - how many times you see each count
    hist[v] += 1       # if there are 3 positions that had 2 counts, histogram at [2] will gain 3 counts

# output
print(f'Size: %d' % (gsize))
print(f'X: %.1f' % (x))
print(f'BACs: %d' % (bacs))
print(f'Min: %d' % (min))
print(f'Max: %d' % (max))
print(f'Counts:')
for i in range(len(hist)): # for each number in the histogram, prints the number and number of counts at that position
    print(i, hist[i])


"""
Size: 1000
X: 5.0
BACs: 5000
Min: 0
Max: 13
Counts:
0 5
1 39
2 88
3 144
4 175
5 150
6 151
7 116
8 59
9 40
10 20
11 5
12 6
13 2
"""
