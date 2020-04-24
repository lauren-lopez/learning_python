#!/usr/bin/env python3

from math import sqrt
import fileinput

# Write a program that computes typical stats
# Count, Min, Max, Mean, Std. Dev, Median
# No, you cannot import any other modules!

numbers = []
count = 0

for line in fileinput.input():
    if line.startswith('#'): continue
    line = line.rstrip()
    numbers.append(float(line))
    count += 1
    numbers.sort()

min = numbers[0]
max = numbers[-1]

sum = 0
sqsum = 0
for n in numbers:
    sum += n
    mean = sum/count
    
# calculate standard deviation
for n in numbers:
    diff = n-mean
    sq = diff**2
    sqsum += sq
    var = sqsum/count
    sd = sqrt(var)
    
# calculate median
if count % 2 == 0:
    med1 = numbers[count//2]
    med2 = numbers[count//2-1]
    median = (med1 + med2)/2

   
print(f'Count: %d' % (count))
print(f'Minimum: %.1f' % (min))
print(f'Maximum: %.1f' % (max))
print(f'Mean: %.3f' % (sum/count))
print(f'Std dev: %.3f' % (sd))
print(f'Median: %.3f' % (median))

"""
python3 stats.py numbers.txt
Count: 10
Minimum: -1.0
Maximum: 256.0
Mean: 29.147789999999997
Std. dev: 75.777
Median 2.35914
"""
