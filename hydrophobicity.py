#!/usr/bin/env python3

# Write a program that computes hydrophobicity in a window
# Let the user choose the method (see below)
# https://en.wikipedia.org/wiki/Hydrophilicity_plot
# https://en.wikipedia.org/wiki/Hydrophobicity_scales

"""
python3 hydrophobicity.py --input proteins.fasta.gz --window 11 --method kd
"""

import argparse
import biotools as bt

# setup
parser = argparse.ArgumentParser(
	description='Computes hydrophobicity within a window')

# required arguments
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='protein file containing sequences')

# optional arguments with default parameters
parser.add_argument('--win', required=False, type=int, default=11,
	metavar='<int>', help='signal peptide window length [%(default)i]')
parser.add_argument('--method', required=False, type=str, default='kd_scale',
	metavar='<str>', help='computation method (options: kd_scale, is_scale, os_scale [%(default)s]')

# finalization
arg = parser.parse_args()

scale = dict()          # This is a dictionary containing all scale dictionaries
kd_scale = {
    'I' : 4.5,  	'V' : 4.2,      'L' : 3.8,	    'F' : 2.8,
    'C' : 2.5,	    'M' : 1.9,	    'A' : 1.8,	    'G' : -0.4,
    'T' : -0.7,	    'S' : -0.8,	    'W' : -0.9,	    'Y' : -1.3,
    'P' : -1.6,	    'H' : -3.2, 	'E' : -3.5,	    'Q' : -3.5,
    'D' : -3.5,	    'N' : -3.5,	    'K' : 3.9,	    'R' : 4.5
    }

is_scale = {
    'I' : -0.31,	'V' : 0.07,     'L' : -0.56,	'F' : -1.13,
    'C' : -0.24,	'M' : -0.23,	'A' : 0.17,	    'G' : 0.01,
    'T' : 0.14,	    'S' : 0.13,	    'W' : -1.85,	'Y' : -0.94,
    'P' : 0.45,	    'H' : 0.17, 	'E' : 2.02,	    'Q' : 0.58,
    'D' : 1.23,	    'N' : 0.42, 	'K' : 0.99,	    'R' : 0.81
    }

os_scale = {
    'I' : -1.12,	'V' : -0.46,    'L' : -1.25,	'F' : -1.71,
    'C' : -0.02,	'M' : -0.67,	'A' : 0.50,	    'G' : 0.11,
    'T' : 0.25,	    'S' : 0.46,	    'W' : -2.09,	'Y' : -0.71,
    'P' : 0.14,	    'H' : 0.11, 	'E' : 3.63,	    'Q' : 0.77,
    'D' : 3.64,	    'N' : 0.85, 	'K' : 2.80,	    'R' : 1.81
    }

isos_scale = {
    'I' : -0.81,	'V' : -0.53,    'L' : -0.69,	'F' : -0.58,
    'C' : 0.22,	    'M' : -0.44,	'A' : 0.33,	    'G' : 0.12,
    'T' : 0.11,	    'S' : 0.33,	    'W' : -0.24,	'Y' : 0.23,
    'P' : -0.31,    'H' : -0.06, 	'E' : 1.61,	    'Q' : 0.19,
    'D' : 2.41,	    'N' : 0.43, 	'K' : 1.81,	    'R' : 1.00
    }

scale['kd_scale'] = dict()          # This puts each dictionary into the larger dictionary
scale['kd_scale'] = kd_scale
scale['is_scale'] = dict()
scale['is_scale'] = is_scale
scale['os_scale'] = dict()
scale['os_scale'] = os_scale
scale['isos_scale'] = dict()
scale['isos_scale'] = isos_scale

win = arg.win

def hydrophobicity(seq, scale_type):        # Calculate hydrophobicity of aa sequence
    hydro_score = 0
    for aa in seq:
        hydro_score += scale[scale_type][aa]
    return hydro_score/len(seq)


for name, seq in bt.read_fasta(arg.file):
    for i in range(len(seq) -win):
        sseq = seq[i:i+win]
        hydro_score = hydrophobicity(sseq, arg.method)
        print(f'{name} {sseq} {hydro_score:.3f}')



