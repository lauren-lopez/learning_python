#!/usr/bin/env python3

import argparse
import biotools as bt

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

# Use a dictionary for the K-D values, and store in biotools
# Use argparse for command line

# setup
parser = argparse.ArgumentParser(
	description='Predicts transmembrane proteins.')

# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='protein file containing sequences')

# optional arguments with default parameters
parser.add_argument('--win1', required=False, type=str, default=8,
	metavar='<int>', help='signal peptide window length [%(default)i]')
parser.add_argument('--win2', required=False, type=int, default=11,
	metavar='<int>', help='hydrophobic peptide window length [%(default)i]')
parser.add_argument('--kd1', required=False, type=float, default=2.5,
	metavar='<float>', help='mean KD for signal peptide [%(default)f]')
parser.add_argument('--kd2', required=False, type=float, default=2.0,
	metavar='<float>', help='mean KD for hydrophobic peptide [%(default)f]')

# finalization
arg = parser.parse_args()


for name, seq in bt.read_fasta(arg.file):
    sseq1 = seq[:30]                                        # first 30 aa for signal peptide
    sseq2 = seq[30:]                                        # everything after the first 30 peptides for hydrophobic region
    if bt.transmembrane(sseq1, arg.win1, arg.kd1) and \
        bt.transmembrane(sseq2, arg.win2, arg.kd2):         # if the condition is met for both regions,
        print(name)



"""
python3 transmembrane.py --file proteins.fasta.gz
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
