#!/usr/bin/env python3

# Modify entropy_fast() however you like to make it faster
# Ideally, your method is faster at all ranges of window size

import math
import time
import random

def entropy_slow(seq, w, th):
	t0 = time.perf_counter()
	low_H_count = 0

	for i in range(len(seq) - w + 1):
		win = seq[i:i+w]
		a, c, g, t = 0, 0, 0, 0
		for nt in win:
			if   nt == 'A': a += 1
			elif nt == 'C': c += 1
			elif nt == 'G': g += 1
			elif nt == 'T': t += 1
		total = a + c + g + t
		h = 0
		pa, pc, pg, pt = a/total, c/total, g/total, t/total

		if a != 0: h -= pa * math.log2(pa)
		if c != 0: h -= pc * math.log2(pc)
		if g != 0: h -= pg * math.log2(pg)
		if t != 0: h -= pt * math.log2(pt)

		if h < th: low_H_count += 1

	t1 = time.perf_counter()
	return low_H_count, t1-t0

def entropy_fast(seq, w, th):
	t0 = time.perf_counter()

	if w < 10:
		low_H_count = 0
		h_win = {}
		for i in range(len(seq) - w + 1):
			win = seq[i:i+w]
			if win not in h_win:
				ha,hg,hc,ht = 0,0,0,0
				pa = win.count("A")/w
				pg = win.count("G")/w
				pc = win.count("C")/w
				pt = win.count("T")/w
				if pa != 0: ha = (pa)*math.log2(pa)
				if pg != 0: hg = (pg)*math.log2(pg)
				if pt != 0: hc = (pt)*math.log2(pt)
				if pc != 0: ht = (pc)*math.log2(pc)
				h = -(ha + hg + hc + ht)
				h_win[win] = h

			if h_win[win] < th: low_H_count += 1

		t1 = time.perf_counter()
		return low_H_count, t1-t0
	elif w < 51:
		low_H_count = 0
		h_win = {}
		for i in range(len(seq) - w + 1):
			win = seq[i:i+w]
			a = win.count("A")
			g = win.count("G")
			c = win.count("C")
			t = win.count("T")

			nt_count = (a,g,c,t)
			if nt_count not in h_win:
				ha,hg,hc,ht = 0,0,0,0
				if a != 0: ha = (a/w)*math.log2(a/w)
				if g != 0: hg = (g/w)*math.log2(g/w)
				if t != 0: hc = (t/w)*math.log2(t/w)
				if c != 0: ht = (c/w)*math.log2(c/w)
				h = -(ha + hg + hc + ht)
				h_win[nt_count] = h

			if h_win[nt_count] < th: low_H_count += 1

		t1 = time.perf_counter()
		return low_H_count, t1-t0
	else:
		t0 = time.perf_counter()
		low_H_count = 0

		a = seq[0:w].count("A")
		g = seq[0:w].count("G")
		c = seq[0:w].count("C")
		t = seq[0:w].count("T")

		ha,hg,hc,ht=0,0,0,0
		if a > 0: ha = (a/w)*math.log2(a/w)
		if g > 0: hg = (g/w)*math.log2(g/w)
		if t > 0: ht = (t/w)*math.log2(t/w)
		if c > 0: hc = (c/w)*math.log2(c/w)
		h = -(ha + hg + hc + ht)
		if h < th : low_H_count += 1

		for i in range(1,len(seq) - w + 1):
			if   seq[i-1]== "A": a -=1
			elif seq[i-1]== "G": g -=1
			elif seq[i-1]== "C": c -=1
			elif seq[i-1]== "T": t -=1

			if   seq[i+w-1]=="A": a +=1
			elif seq[i+w-1]=="G": g +=1
			elif seq[i+w-1]=="C": c +=1
			elif seq[i+w-1]=="T": t +=1

			ha,hg,hc,ht=0,0,0,0

			if a > 0: ha = (a/w)*math.log2(a/w)
			if g > 0: hg = (g/w)*math.log2(g/w)
			if t > 0: ht = (t/w)*math.log2(t/w)
			if c > 0: hc = (c/w)*math.log2(c/w)
			h = -(ha + hg + hc + ht)
			if h < th : low_H_count += 1

		t1 = time.perf_counter()
		return low_H_count, t1-t0

# create a random chromosome
seq = []
alph = ['A', 'C', 'G', 'T']
for i in range(int(1e5)):
	seq.append(alph[random.randint(0,3)])
seq = ''.join(seq)

# test speed at various word sizes
W = [2, 7, 15, 100, 1000]
for w in W:
	cs, ts = entropy_slow(seq, w, 1)
	cf, tf = entropy_fast(seq, w, 1)
	assert(cs == cf)
	print(ts / tf)
