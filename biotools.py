#!/usr/bin/env python3

import sys
import gzip
import random

def read_fasta(filename):                   # Reads FASTA files
	name = None
	seqs = []
	
	fp = None
	if filename == '-':
		fp = sys.stdin
	elif filename.endswith('.gz'):
		fp = gzip.open(filename, 'rt')
	else:
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()


def gc(seq):                                # Calculates the average GC content of a sequence
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)


def gc_skew(seq):                           # Calculates the GC Skew:(g-c)/(g+c)
    g, c = 0, 0
    for nt in seq:
        if nt   == 'G': g += 1
        elif nt == 'C': c += 1
    return((g-c)/(g+c))


def randseq(length, gc):                    # Generates random DNA sequences for a defined GC content
    seq = []
    for nt in range(length):
        r = random.random()
        if r < gc: 
            r = random.random()
            if r < 0.5: seq.append('G')
            else:       seq.append('C')
        else:
            r = random.random()
            if r < 0.5: seq.append('A')
            else:       seq.append('T')
    return ''.join(seq)

"""
def kd(seq):                                # Calculates the mean KD score for an amino acid sequence
    score = 0       # count up individual hydropathy scores of all amino acids
    for aa in seq:
        if   aa == 'I': score += 4.5
        elif aa == 'V': score += 4.2
        elif aa == 'L': score += 3.8
        elif aa == 'F': score += 2.8
        elif aa == 'C': score += 2.5
        elif aa == 'M': score += 1.9
        elif aa == 'A': score += 1.8
        elif aa == 'G': score -= 0.4
        elif aa == 'T': score -= 0.7
        elif aa == 'S': score -= 0.8
        elif aa == 'W': score -= 0.9
        elif aa == 'Y': score -= 1.3
        elif aa == 'P': score -= 1.6
        elif aa == 'H': score -= 3.2
        elif aa == 'E': score -= 3.5
        elif aa == 'Q': score -= 3.5
        elif aa == 'D': score -= 3.5
        elif aa == 'N': score -= 3.5
        elif aa == 'K': score -= 3.9
        elif aa == 'R': score -= 4.5
    return (score/len(seq))
"""

def kd(seq):                                # Calculates the mean KD score for an amino acid sequence
    score = 0       # count up individual hydropathy scores of all amino acids
    kd_scale = {
        'I' : '4.5',	'V' : '4.2',    'L' : '3.8',	'F' : '2.8',
        'C' : '2.5',	'M' : '1.9',	'A' : '1.8',	'G' : '-0.4',
        'T' : '-0.7',	'S' : '-0.8',	'W' : '-0.9',	'Y' : '-1.3',
        'P' : '-1.6',	'H' : '-3.2',	'E' : '-3.5',	'Q' : '-3.5',
        'D' : '-3.5',	'N' : '-3.5',	'K' : '3.9',	'R' : '4.5',
    }
    
    for aa in seq:
        if aa in kd_scale: score += float(kd_scale[aa])   
    return (score/len(seq))


def transmembrane(seq, win, threshold):     # Identifies proteins that are predicted to be transmembrane
    for i in range(len(seq)-win+1):                         # aa sequence, window size, KD threshold
        sseq = seq[i : i+win]                               # sseq is the aa within the window
        if kd(sseq) > threshold and 'P' not in sseq:
            return  True                                    # if a True is found, the function exits here
    return  False                                           # if a True is not found, the function returns False


def longest_orf(seq):                   # Finds the longest ORF in an mRNA sequence and returns the sequence
# find start codon (any reading frame)
    atgs = []
    for i in range(len(seq) -2):
        codon = seq[i:i+3]
        if codon == 'ATG': atgs.append(i)

# find in-frame stop codon    
    max_len = 0
    max_seq = None
    for start in atgs:
        stop = None
        for i in range(start, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon == 'TAA' or codon == 'TGA' or codon == 'TAG':
                stop = i
                break
        if stop != None:
            cds_len = stop - start +3
            if cds_len > max_len:
                max_len = cds_len
                max_seq = seq[start : start+cds_len]
    if max_seq == None: return None
    return max_seq

"""
def translate(seq):                         # Translates mRNA sequences to aa
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
"""  

def translate(seq):

    gcode = {
        'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
        'ACA' : 'T',	'ACC' : 'T',	'ACG' : 'T',	'ACT' : 'T',
        'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',	'AGT' : 'S',
        'ATA' : 'I',	'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',
        'CAA' : 'Q',	'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',
        'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',	'CCT' : 'P',
        'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
        'CTA' : 'L',	'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',
        'GAA' : 'E',	'GAC' : 'D',	'GAG' : 'E',	'GAT' : 'D',
        'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
        'GGA' : 'G',	'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',
        'GTA' : 'V',	'GTC' : 'V',	'GTG' : 'V',	'GTT' : 'V',
        'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',	'TAT' : 'Y',
        'TCA' : 'S',	'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',
        'TGA' : '*',	'TGC' : 'C',	'TGG' : 'W',	'TGT' : 'C',
        'TTA' : 'L',	'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',
    }

    pro_seq = []
    assert len(seq) % 3 == 0
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if codon in gcode: pro_seq.append(gcode[codon])
        else:              pro_seq.append('X')

    return ''.join(pro_seq)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    