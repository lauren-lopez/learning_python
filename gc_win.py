#!/usr/bin/env python3

# Write a program that computes the GC fraction of a DNA sequence in a window
# Window size is 11 nt
# Output with 4 significant figures using whichever method you prefer

seq = 'ACGACGCAGGAGGAGAGTTTCAGAGATCACGAATACATCCATATTACCCAGAGAGAG'
w = 11
s = 0 #start position
gc_count = 0

for i in range(s, len(seq)-w+1): #doesn't work - prints total number GC
    for c in seq[s:s+w]:
        if c == 'G' or c == 'C': gc_count += 1
        gc_frac = gc_count/w
    print('%d %s %.4f' % (i, seq[s:s+w], gc_frac))
    s += 1
    gc_count = 0


#Failed attempts below - kept for future reference 
#for i in range(s, len(seq)-w+1): #doesn't work - prints total number GC
#    for c in seq:
#        if c == 'G' or c == 'C': gc_count += 1
#    print(i, seq[s:s+w], gc_count)
#    s += 1
#    gc_count = 0

#for i in range(s, len(seq)-w+1): # doesn't work - gc count always 0
#    if seq[i:(i+w)] == 'G' or seq[i:(i+w)] == 'C': gc_count += 1 
#    print(i, seq[s:s+w], gc_count)
#    s += 1
#    gc_count = 0

#for i in range(s, len(seq)-w+1): # doesn't work - gc count always 0
#    if seq[i:i+w] == 'G' or seq[i:i+w] == 'C': gc_count += 1 
#    print(i, seq[s:s+w], gc_count)
#    s += 1
#    gc_count = 0

#for i in range(s, len(seq)-w+1): #doesn't work - only checks first letter in window
#    if seq[i] == 'G' or seq[i] == 'C': gc_count += 1 
#    print(i, seq[s:s+w], gc_count)
#    s += 1
#    gc_count = 0 #resets GC count within 11-nt window back to 0 after each loop
    
#for i in range(s, len(seq)-w+1): #works, but pretty sure method is wrong
#    if seq[i]    == 'G' or seq[i]    == 'C': gc_count += 1
#    if seq[i+1]  == 'G' or seq[i+1]  == 'C': gc_count += 1
#    if seq[i+2]  == 'G' or seq[i+2]  == 'C': gc_count += 1
#    if seq[i+3]  == 'G' or seq[i+3]  == 'C': gc_count += 1
#    if seq[i+4]  == 'G' or seq[i+4]  == 'C': gc_count += 1
#    if seq[i+5]  == 'G' or seq[i+5]  == 'C': gc_count += 1
#    if seq[i+6]  == 'G' or seq[i+6]  == 'C': gc_count += 1
#    if seq[i+7]  == 'G' or seq[i+7]  == 'C': gc_count += 1
#    if seq[i+8]  == 'G' or seq[i+8]  == 'C': gc_count += 1
#    if seq[i+9]  == 'G' or seq[i+9]  == 'C': gc_count += 1
#    if seq[i+10] == 'G' or seq[i+10] == 'C': gc_count += 1
#    gc_frac = gc_count/w
#    print('%d %s %.4f' % (i, seq[s:s+w], gc_frac))
#    s += 1
#    gc_count = 0 #resets GC count within 11-nt window back to 0 after each loop

"""
0 ACGACGCAGGA 0.6364
1 CGACGCAGGAG 0.7273
2 GACGCAGGAGG 0.7273
3 ACGCAGGAGGA 0.6364
4 CGCAGGAGGAG 0.7273
5 GCAGGAGGAGA 0.6364
6 CAGGAGGAGAG 0.6364
7 AGGAGGAGAGT 0.5455
8 GGAGGAGAGTT 0.5455
9 GAGGAGAGTTT 0.4545
10 AGGAGAGTTTC 0.4545
11 GGAGAGTTTCA 0.4545
12 GAGAGTTTCAG 0.4545
13 AGAGTTTCAGA 0.3636
14 GAGTTTCAGAG 0.4545
15 AGTTTCAGAGA 0.3636
16 GTTTCAGAGAT 0.3636
17 TTTCAGAGATC 0.3636
18 TTCAGAGATCA 0.3636
19 TCAGAGATCAC 0.4545
20 CAGAGATCACG 0.5455
21 AGAGATCACGA 0.4545
22 GAGATCACGAA 0.4545
23 AGATCACGAAT 0.3636
24 GATCACGAATA 0.3636
25 ATCACGAATAC 0.3636
26 TCACGAATACA 0.3636
27 CACGAATACAT 0.3636
28 ACGAATACATC 0.3636
29 CGAATACATCC 0.4545
30 GAATACATCCA 0.3636
31 AATACATCCAT 0.2727
32 ATACATCCATA 0.2727
33 TACATCCATAT 0.2727
34 ACATCCATATT 0.2727
35 CATCCATATTA 0.2727
36 ATCCATATTAC 0.2727
37 TCCATATTACC 0.3636
38 CCATATTACCC 0.4545
39 CATATTACCCA 0.3636
40 ATATTACCCAG 0.3636
41 TATTACCCAGA 0.3636
42 ATTACCCAGAG 0.4545
43 TTACCCAGAGA 0.4545
44 TACCCAGAGAG 0.5455
45 ACCCAGAGAGA 0.5455
46 CCCAGAGAGAG 0.6364
"""
