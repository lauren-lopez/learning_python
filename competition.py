#!/usr/bin/env python3

# Modify entropy_fast() however you like to make it faster
# Ideally, your method is faster at all ranges of window size

import math
import time
import random

def entropy_slow(seq, w, th):
    t0 = time.perf_counter()
    low_H_count = 0

    for i in range(len(seq) - w + 1):
        win = seq[i:i+w]
        a, c, g, t = 0, 0, 0, 0
        for nt in win:                                          # count up nts in a window as variables
            if   nt == 'A': a += 1
            elif nt == 'C': c += 1
            elif nt == 'G': g += 1
            elif nt == 'T': t += 1
        total = a + c + g + t
        h = 0
        pa, pc, pg, pt = a/total, c/total, g/total, t/total     # proportion of each nt in window

        if a != 0: h -= pa * math.log2(pa)
        if c != 0: h -= pc * math.log2(pc)
        if g != 0: h -= pg * math.log2(pg)
        if t != 0: h -= pt * math.log2(pt)

        if h < th: low_H_count += 1                             # counts up the number of windows that have low entropy (<h)

    t1 = time.perf_counter()
    return low_H_count, t1-t0

def entropy_fast(seq, w, th):
    t0 = time.perf_counter()
    low_H_count = 0
    win_entropy = {}                                # creates a dictionary containing window sequence and entropy value
    
    for i in range(len(seq) - w + 1):
        win = seq[i:i+w]
        if win in win_entropy: 
            h = win_entropy[win]
            if h < th: 
                low_H_count += 1
                continue
        else:
            count = {}                                  # nt proportions go into dictionary instead of being counted up in variables
            for nt in win:
                if nt in count: count[nt] += 1/w        # 1/w skips a few steps to add the nt proportions into the dictionary
                else:           count[nt] = 1/w         # instead of counting up nt then dividing by totals

            h = 0
            for p in count.values():
                if p != 0: h -= p * math.log2(p)
                
            if h < th: low_H_count += 1
            win_entropy[win] = h

    t1 = time.perf_counter()
    return low_H_count, t1-t0
    

# create a random chromosome
seq = []
alph = ['A', 'C', 'G', 'T']
for i in range(int(1e5)):
    seq.append(alph[random.randint(0,3)])
seq = ''.join(seq)

# test speed at various window sizes
W = [3, 4, 7, 10, 20, 100]
for w in W:
    cs, ts = entropy_slow(seq, w, 1)            
    cf, tf = entropy_fast(seq, w, 1)             
    assert(cs == cf)                        # ensures both methods counted the same number of low entropy windows
    print(ts/tf, w)                         # prints how much faster the 'fast' method is compared to the 'slow'

"""
def entropy_fast(seq, w, th):
    t0 = time.perf_counter()
    low_H_count = 0
    win_entropy = {}                                # creates a dictionary containing window sequence and entropy value
    
    for i in range(len(seq) - w + 1):
        win = seq[i:i+w]
        if win in win_entropy: 
            h = win_entropy[win]
            if h < th: 
                low_H_count += 1
                continue
        else:
            count = {}                                  # nt proportions go into dictionary instead of being counted up in variables
            for nt in win:
                if nt in count: count[nt] += 1/w        # 1/w skips a few steps to add the nt proportions into the dictionary
                else:           count[nt] = 1/w         # instead of counting up nt then dividing by totals

            h = 0
            for p in count.values():
                if p != 0: h -= p * math.log2(p)
                
            if h < th: low_H_count += 1
            win_entropy[win] = h

    t1 = time.perf_counter()
    return low_H_count, t1-t0
"""
