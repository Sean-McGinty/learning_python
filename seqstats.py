#!/usr/bin/env python3

import fileinput

# Write a program that computes typical sequence stats
# No, you cannot import any other modules!
# Use rand_seq to generate the sequences
# Expected output is shown below

seq_size = []
l_count = 0
total_seq = ""
a_count = 0
c_count = 0
g_count = 0
t_count = 0
seq_count = 0
for line in fileinput.input():
    if line.startswith(">"): continue
    seq=line.rstrip()
    l_count+=len(seq)
    for nt in seq:
        seq_count+=1
        if nt=="A":
            a_count += 1
        elif nt=="G":
            g_count += 1
        elif nt=="T":
            t_count += 1
        else:
            c_count += 1
    seq_size.append(len(seq))
seq_size.sort()
sum = 0
i = len(seq_size)
while sum <= l_count/2:
    i -= 1
    sum += seq_size[i]
print(f'Number of Sequences: {len(seq_size)}')
print(f'Number of Letters: {l_count}')
print(f'Minimum Length: {seq_size[0]}')
print(f'Maximum Length: {seq_size[-1]}')
print(f'N50: {seq_size[i]}')
print(f'A={a_count/l_count:.4f} C={c_count/l_count:.4f} G={g_count/l_count:.4f} T={t_count/l_count:.4f}')
"""
python3 rand_seq.py 100 100 100000 0.35 | python3 seqstats.py
Number of sequences: 100
Number of letters: 4957689
Minimum length: 219
Maximum length: 99853
N50: 67081
Composition: A=0.325 C=0.175 G=0.175 T=0.325
"""
