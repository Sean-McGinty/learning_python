#!/usr/bin/env python3

# Write a program that masks areas of low complexity sequence
# Use argparse for command line arguments (see example below)
# Use read_fasta() from biotools.py
import gzip
import sys
import biotools as bt
import argparse

parser = argparse.ArgumentParser(
	description='Low complexity sequence masker.')

# required arguments
parser.add_argument('--input', required=True, type=str,
	metavar='<path>', help='fasta file')

# optional arguments with default parameters
parser.add_argument('--window', required=False, type=int, default=15,
	metavar='<int>', help='ptional integer argument [%(default)i]')
parser.add_argument('--threshold', required=False, type=float,default=1.100000,
	metavar='<float>', help='entropy threshold [%(default)f]')
parser.add_argument('--lowercase', action='store_true',
	help='report lower case instead of N')
arg = parser.parse_args()

#if arg.lowercase:letter = "n"
#else:            letter = "N"

for name, seq in bt.read_fasta(arg.input):
	print(f'>{name}')
	sequence=[]
	for nt in seq: sequence.append(nt)
	for i in range(0,len(seq)-arg.window+1,1):
		win=seq[i:i+arg.window]
		if bt.entropy(win) < arg.threshold:
			if arg.lowercase:
				sequence[i:i+arg.window]=win.lower()
			else: sequence[i:i+arg.window]= "N"*arg.window
	print("".join(sequence))

"""
python3 entropy_filter.py --help
usage: entropy_filter.py [-h] --input <path> [--window <int>]
                         [--threshold <float>] [--lowercase]

Low complexity sequence masker.

optional arguments:
  -h, --help           show this help message and exit
  --input <path>       fasta file
  --window <int>       optional integer argument [15]
  --threshold <float>  entropy threshold [1.100000]
  --lowercase          report lower case instead of N


python3 entropy_filter.py --input genome.fa.gz | head -20
>I
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA
GCCTAAGCCTAAAAAATTGAGATAAGAAAACATTTTACTTTTTCAAAATTGTTTTCATGC
TAAATTCAAAACNNNNNNNNNNNNNNNAAGCTTCTAGATATTTGGCGGGTACCTCTAATT
TTGCCTGCCTGCCAACCTATATGCTCCTGTGTTTAGGCCTAATACTAAGCCTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAA
GCCTAAGACTAAGCCTAAGACTAAGCCTAATACTAAGCCTAAGCCTAAGACTAAGCCTAA
GCCTAAAAGAATATGGTAGCTACAGAAACGGTAGTACACTCTTCTGNNNNNNNNNNNNNN
NTGCAATTTTTATAGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAAT
TGCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAAACANNNNNNNNNNNNNNNGAAAT
GAATATCGTAGCTACAGAAACGGTTGTGCACTCATCTGAAANNNNNNNNNNNNNNNNNNN
NNGCACTTTGTGCAGAATTCTTGATTCTTGATTCTTGCAGAAATTTGCAAGAAAATTCGC
"""
