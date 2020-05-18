#!/usr/bin/env python3

import gzip
import sys
import biotools as bt
import argparse

# Use argparse
# Write a program that translates an mRNA
# Assume the protein encoded is the longest ORF

def longest_orf(seq):
    # find all ATGs (start codon)
    assert len(seq) > 0
    atgs = []
    for i in range(len(seq) -2):
        if seq[i:i+3] == 'ATG': atgs.append(i)              # position number of each START codon
    
    # for each ATG, find nearest in-frame STOP
    # check if longest
    max_len = 0
    max_seq = None
    for atg in atgs:
        stop = None
        for i in range(atg, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                stop = i                                    # position number of the first STOP codon after each START
                break                                       # Only need to find the first, not all
        if stop != None:
            cds_len = stop - atg +3                         # cds = coding region *Make sure ends with a STOP CODON!*
            if cds_len > max_len:                           # checks if sequence is longest
                max_len = cds_len
                max_seq = seq[atg:atg+cds_len]              # max_seq is the seq with the longest orf length
    
    # translate longest ORF into protein
    if max_seq == None: return None
    return translate(max_seq)
    
def translate(seq):
    assert(len(seq) % 3 == 0)                               # ensures the seq len is divisible by 3, i.e. full codons
    pro = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if   codon == 'AAA': pro.append('K')
        elif codon == 'AAC': pro.append('N')
        elif codon == 'AAG': pro.append('K')
        elif codon == 'AAT': pro.append('N')
        elif codon == 'ACA': pro.append('T')
        elif codon == 'ACC': pro.append('T')
        elif codon == 'ACG': pro.append('T')
        elif codon == 'ACT': pro.append('T')
        elif codon == 'AGA': pro.append('R')
        elif codon == 'AGC': pro.append('S')
        elif codon == 'AGG': pro.append('R')
        elif codon == 'AGT': pro.append('S')
        elif codon == 'ATA': pro.append('I')
        elif codon == 'ATC': pro.append('I')
        elif codon == 'ATG': pro.append('M')
        elif codon == 'ATT': pro.append('I')
        
        elif codon == 'CAA': pro.append('Q')
        elif codon == 'CAC': pro.append('H')
        elif codon == 'CAG': pro.append('Q')
        elif codon == 'CAT': pro.append('H')
        elif codon == 'CCA': pro.append('P')
        elif codon == 'CCC': pro.append('P')
        elif codon == 'CCG': pro.append('P')
        elif codon == 'CCT': pro.append('P')
        elif codon == 'CGA': pro.append('R')
        elif codon == 'CGC': pro.append('R')
        elif codon == 'CGG': pro.append('R')
        elif codon == 'CGT': pro.append('R')
        elif codon == 'CTA': pro.append('L')
        elif codon == 'CTC': pro.append('L')
        elif codon == 'CTG': pro.append('L')
        elif codon == 'CTT': pro.append('L')
        
        elif codon == 'GAA': pro.append('E')
        elif codon == 'GAC': pro.append('D')
        elif codon == 'GAG': pro.append('E')
        elif codon == 'GAT': pro.append('D')
        elif codon == 'GCA': pro.append('A')
        elif codon == 'GCC': pro.append('A')
        elif codon == 'GCG': pro.append('A')
        elif codon == 'GCT': pro.append('A')
        elif codon == 'GGA': pro.append('G')
        elif codon == 'GGC': pro.append('G')
        elif codon == 'GGT': pro.append('G')
        elif codon == 'GGG': pro.append('G')
        elif codon == 'GTA': pro.append('V')
        elif codon == 'GTC': pro.append('V')
        elif codon == 'GTG': pro.append('V')
        elif codon == 'GTT': pro.append('V')
        
        elif codon == 'TAA': pro.append('*')
        elif codon == 'TAC': pro.append('Y')
        elif codon == 'TAG': pro.append('*')
        elif codon == 'TAT': pro.append('Y')
        elif codon == 'TCA': pro.append('S')
        elif codon == 'TCC': pro.append('S')
        elif codon == 'TCT': pro.append('S')
        elif codon == 'TCT': pro.append('S')
        elif codon == 'TGA': pro.append('*')
        elif codon == 'TGC': pro.append('C')
        elif codon == 'TGG': pro.append('W')
        elif codon == 'TGT': pro.append('C')
        elif codon == 'TTA': pro.append('L')
        elif codon == 'TTC': pro.append('F')
        elif codon == 'TTG': pro.append('L')
        elif codon == 'TTT': pro.append('F')
        else:                pro.append('X')
    return ''.join(pro)

    
for name, seq in bt.read_fasta('mRNA.fa.gz'):
    pro = longest_orf(seq)
    if pro != None:
        print(f'>{name}')
        print(pro)



"""
python3 translate_mRNA.py --file ../Lesson05/transcripts.fasta.gz
>CBG00001.1
MTFCENKNLPKPPSDRCQVVVISILSMILDFYLKYNPDKHWAHLFYGASPILEILVIFGMLANSVYGNKLAMFACVLDLVSGVFCLLTLPVISVAENATGVRLHLPYISTFHSQFSFQVSTPVDLFYVATFLGFVSTILILLFLILDALKFMKLRKLRNEDLEKEKKMNPIEKV*
>CBG00006.1
MNGVEKVNKYFDIKDKRDFLYHFGFGVDTLDIKAVFGDTKFVCTGGSPGRFKLYAEWFAKETSIPCSENLSRSDRFVIYKTGPVCWINHGMGTPSLSIMLVESFKLMHHAGVKNPTFIRLGTSGGVGVPPGTVVVSTGAMNAELGDTYVQVIAGKRIERPTQLDATLREALCAVGKEKNIPVETGKTMCADDFYEGQMRLDGYFCDYEEEDKYAFLRKLNSLGVRNIEMESTCFASFTCRAGFPSAIVCVTLLNRMDGDQVQIDKEKYIEYEERPFRLVTAYIRQQTGV*
etc.
"""
