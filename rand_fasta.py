#!/usr/bin/env python3

import gzip
import sys
import math
import random

# Write a program that finds creates random fasta files
# Create a function that makes random DNA sequences
# Parameters include length and frequencies for A, C, G, T
# Command line: python3 rand_fasta.py <count> <min> <max> <a> <c> <g> <t>


def rand_seq(lenseq, a, g, t, c):
    dna_seq = []
    for i in range(lenseq):
        r=random.random()
        if   r<a:       dna_seq.append("A")
        elif r<a+c:     dna_seq.append("C")
        elif r<a+c+t:   dna_seq.append("T")
        elif r<a+c+t+g: dna_seq.append("G")
    return''.join(dna_seq)

count = int(sys.argv[1])
min =   int(sys.argv[2])
max =   int(sys.argv[3])
a =     float(sys.argv[4])
c =     float(sys.argv[5])
g =     float(sys.argv[6])
t =     float(sys.argv[7])

for i in range(count):
    x = random.randint(min,max)
    dna = rand_seq(x, a, g, t, c)
    print(f'>{i}')
    print(dna)
