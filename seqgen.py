#!/usr/bin/env python3

import argparse
import random
import biotools as bt

# setup
parser = argparse.ArgumentParser(
	description='Generates random FASTA files.')
# required arguments
parser.add_argument('--count', required=True, type=int,
	metavar='<int>', help='number of sequences')
parser.add_argument('--min', required=True, type=int,
	metavar='<int>', help='minimum sequence length')
parser.add_argument('--max', required=True, type=int,
	metavar='<int>', help='maximum sequence length')
parser.add_argument('--gc', required=True, type=float,
	metavar='<float>', help='gc composition')
parser.add_argument('--prefix', required=False, type=str, default='seq',
	metavar='<str>', help='sequence name prefix')

# switches
parser.add_argument('--verbose', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()

for i in range(arg.count):
    print(f'>{arg.prefix}-{i}')
    slen = random.randint(arg.min, arg.max)
    print(bt.randseq(slen, arg.gc))