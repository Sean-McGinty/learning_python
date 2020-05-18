#!/usr/bin/env python3

import sys
import gzip
import random
import math

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

def gc(seq):
	count = 0
	for nt in seq:
		if nt == 'G' or nt == 'C':
			count += 1
	return count / len(seq)

def randseq(l, gc):
	dna = []
	for i in range(l):
		r = random.random()
		if r < gc:
			r = random.random()
			if r < .5: dna.append("G")
			else:      dna.append("C")
		else:
			r = random.random()
			if r < .5: dna.append("A")
			else:      dna.append("T")
	return "".join(dna)


def kd (seq):
	kd=0
	for aa in seq:
		if   aa == "I": kd+=4.5
		elif aa == "V": kd+=4.2
		elif aa == "L": kd+=3.8
		elif aa == "F": kd+=2.8
		elif aa == "C": kd+=2.5
		elif aa == "M": kd+=1.9
		elif aa == "A": kd+=1.8
		elif aa == "G": kd+=-.4
		elif aa == "T": kd+=-.7
		elif aa == "S": kd+=-.8
		elif aa == "W": kd+=-.9
		elif aa == "Y": kd+=-1.3
		elif aa == "P": kd+=-1.6
		elif aa == "H": kd+=-3.2
		elif aa == "E": kd+=-3.5
		elif aa == "Q": kd+=-3.5
		elif aa == "D": kd+=-3.5
		elif aa == "N": kd+=-3.5
		elif aa == "K": kd+=-3.9
		elif aa == "R": kd+=-4.5
	return kd/len(seq)

def hydrophobic(seq,size,th):
	for i in range(0,len(seq) - size + 1, 1):
		signal=seq[i:i+size]
		if not ("P" in signal) and kd(signal) > th:
			return True
	return False

def skew(seq):
	s=(seq.count("G") - seq.count("C"))/(seq.count("G") + seq.count("C"))
	return s

def rev_comp(seq):
	comp = ""
	for nt in seq:
		if   nt == "A": comp="T" + comp
		elif nt == "G": comp="C" + comp
		elif nt == "C": comp="G" + comp
		elif nt == "T": comp="A" + comp
	return comp

def max_orf(seq):
	atgs=[]
	max_length = 0
	orf = None

	for i in range(len(seq)):
		codon = seq[i:i+3]
		if codon == "ATG":
			atgs.append(i)
	for start_positon in atgs:
		codon_count=0
		for i in range(start_positon,len(seq),3):
			codon=seq[i:i+3]
			codon_count+=1
			if codon == "TAA" or codon == "TAG" or codon == "TGA":
				if codon_count > max_length:
					max_length = codon_count
					orf = seq[start_positon:start_positon + codon_count*3]
				break
	return orf

def DNA_to_AA(orf):
    protein_seq=[]
    for i in range(0,len(orf),3):
        codon = orf[i:i+3]
        if   codon == "TTT" or codon == "TTC":                                                                         protein_seq.append("F")
        elif codon == "TTA" or codon == "TTG" or codon == "CTT" or codon == "CTC" or codon == "CTA" or codon == "CTG": protein_seq.append("L")
        elif codon == "ATT" or codon == "ATC" or codon == "ATC":                                                       protein_seq.append("I")
        elif codon == "ATG":                                                                                           protein_seq.append("M")
        elif codon == "GTT" or codon == "GTC" or codon == "GTA" or codon == "GTG":                                     protein_seq.append("V")
        elif codon == "TCT" or codon == "TCC" or codon == "TCA" or codon == "TCG" or codon == "AGT" or codon == "AGC": protein_seq.append("S")
        elif codon == "CCT" or codon == "CCC" or codon == "CCA" or codon == "CCG":                                     protein_seq.append("P")
        elif codon == "ACT" or codon == "ACC" or codon == "ACA" or codon == "ACG":                                     protein_seq.append("T")
        elif codon == "GCT" or codon == "GCC" or codon == "GCA" or codon == "GCG":                                     protein_seq.append("A")
        elif codon == "TAT" or codon == "TAC":                                                                         protein_seq.append("Y")
        elif codon == "CAT" or codon == "CAC":                                                                         protein_seq.append("H")
        elif codon == "CAA" or codon == "CAG":                                                                         protein_seq.append("Q")
        elif codon == "AAT" or codon == "AAC":                                                                         protein_seq.append("N")
        elif codon == "AAA" or codon == "AAG":                                                                         protein_seq.append("K")
        elif codon == "GAT" or codon == "GAC":                                                                         protein_seq.append("D")
        elif codon == "GAA" or codon == "GAG":                                                                         protein_seq.append("E")
        elif codon == "TGT" or codon == "TGC":                                                                         protein_seq.append("C")
        elif codon == "TGG":                                                                                           protein_seq.append("W")
        elif codon == "CGT" or codon == "CGC" or codon == "CGA" or codon == "CGG" or codon == "AGA" or codon == "AGG": protein_seq.append("R")
        elif codon == "GGT" or codon == "GGC" or codon == "GGA" or codon == "GGG":                                     protein_seq.append("G")
    return ''.join(protein_seq)

def entropy(seq):
	l = len(seq)

	a = seq.count("A")
	g = seq.count("G")
	c = seq.count("C")
	t = seq.count("T")

	ha = 0
	hg = 0
	hc = 0
	ht = 0

	if a != 0: ha = (a/l)*math.log2(a/l)
	if g != 0: hg = (g/l)*math.log2(g/l)
	if t != 0: hc = (t/l)*math.log2(t/l)
	if c != 0: ht = (c/l)*math.log2(c/l)

	h = ha + hg + hc + ht

	return -h
