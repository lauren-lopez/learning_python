#!/usr/bin/env python3

import gzip
import sys

# Write a program that predicts if a protein is trans-membrane
# Trans-membrane proteins have the following properties
#	Signal peptide: https://en.wikipedia.org/wiki/Signal_peptide
#	Hydrophobic regions(s): https://en.wikipedia.org/wiki/Transmembrane_protein
#	No prolines (alpha helix)
# Hydrophobicity is measued via Kyte-Dolittle
#	https://en.wikipedia.org/wiki/Hydrophilicity_plot
# For our purposes:
#	Signal peptide is 8 aa long, KD > 2.5, first 30 aa
#	Hydrophobic region is 11 aa long, KD > 2.0, after 30 aa

def read_fasta(filename):
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


def kd(seq):
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

def transmembrane(seq, win, threshold):                     # aa sequence, window size, KD threshold
    for i in range(len(seq)-win+1):
        sseq = seq[i : i+win]                               # sseq is the aa within the window
        if kd(sseq) > threshold and 'P' not in sseq:
            return  True                                    # if a True is found, the function exits here
    return  False                                           # if a True is not found, the function returns False

for name, seq in read_fasta(sys.argv[1]):
    sseq1 = seq[:30]                                        # first 30 aa for signal peptide
    sseq2 = seq[30:]                                        # everything after the first 30 peptides for hydrophobic region
    if transmembrane(sseq1, 8, 2.5) and transmembrane(sseq2, 11, 2.0):      # if the condition is met for both regions,
        print(name)                                                         # print the name of the protein



# Function works, but not written in the correct way for the assignment
"""
# sig_region = how many aa are in the signal peptide region of the sequence
# win1 = how long the signal peptide should be
# win2 = how long the hydrophobic region should be
# threshold1 = the minimum KD value of the signal peptide
# threshold2 = the minimum KD value of the hydrophobic region
def transmembrane(sig_region, win1, win2, threshold1, threshold2):
    protein_list = []    
    for name, seq in read_fasta(sys.argv[1]):
        sseq1 = seq[:sig_region]                        # sseq1 is the first 30 aa range for the signal peptide
        for i in range(len(sseq1)-win1+1):
            sigpep = sseq1[i : i+win1]                  # sigpep is every 8aa window within the 30aa range
            sigpep_found = False
            if kd(sigpep) > threshold1:                 # sigpep_found is True if at least one 8aa region is found that does not contain a P
                if 'P' not in sigpep:                   # loop breaks after this condition is met - no need to keep searching
                    sigpep_found = True
                    break
        sseq2 = seq[sig_region:]                        # sseq2 is the remainder of the sequence after the signal peptide range
        for j in range(len(sseq2)-win2+1):
            hydro = sseq2[j : j+win2]                   # hydro is every 11aa window within the remainder of the sequence
            hydro_found = False
            if kd(hydro) > threshold2:
                if 'P' not in hydro:
                    hydro_found = True
                    break
        if sigpep_found and hydro_found:                # if a signal peptide and a hydrophobic region are both found, 
            protein_list.append(name)                   # the protein name is added to the list of transmembrane proteins
    return(protein_list)

prot = transmembrane(30, 8, 11, 2.5, 2.0)
print(f'Predicted Transmembrane Proteins: {prot} {len(prot)}')
"""


# First attempt - code did not work
"""
# Find all possible signal peptides and create a list of protein that have a signal peptide
sigpep_len = 8
sigpep_prot = []
for name, seq in read_fasta(sys.argv[1]):
    signal_range = seq[0:30]                    # The signal peptide can be within the first 30 aa
    for i in range(len(signal_range) -sigpep_len+1):
        sigpep = seq[i : i+sigpep_len]          # sigpep is the 8aa window within the first 30aa range
        if kd(sigpep) > 2.5:                    # if the mean kd is >2.5 and the sequence does not contain P, 
            if 'P' not in sigpep:               # then we have found a possible signal peptide
                if name not in sigpep_prot:     # if the protein name is not already in the list, 
                    sigpep_prot.append(name)    # append it to the list of proteins that have a possible signal peptide    
#print(sigpep_prot)

# Within the list of proteins with a signal peptide, find the proteins that have a hydrophobic region
hydro_len = 11
transmembrane = []
for prot in sigpep_prot:
    hydro_range = prot[30:]                      # The hydrophobic region is after the first 30 aa
    for i in range(len(hydro_range) -hydro_len+1):
        hydro = seq[i : i+hydro_len]            # hydro is the 11aa window
        if kd(hydro) > 2.0:
            if 'P' not in hydro:
                if name not in transmembrane:
                    transmembrane.append(name)
print(transmembrane)
"""

"""
python3 transmembrane.py proteins.fasta.gz
18w
Dtg
Krn
Lac
Mcr
PRY
Pxt
Pzl
QC
Ror
S1P
S2P
Spt
apn
bai
bdl
bou
bug
cue
drd
ft
grk
knk
ksh
m
nac
ort
rk
smo
thw
tsg
waw
zye
"""
