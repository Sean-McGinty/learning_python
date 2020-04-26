#!/usr/bin/env python3

# Write a Shannon entropy calculator: H = sum(pi * log(pi))
# Use fileinput to get the data from nucleotides.txt
# Make sure that the values are probabilities
# Make sure that the distribution sums to 1
# Report with 3 decimal figures
import fileinput
import math
data = []
sum=0
for line in fileinput.input():
    if line.startswith('#'): continue # same as above
    vals=line.split()
    print(vals)
    v=float(vals[1])
    print(v)
    assert(v>=0 and v<=1)
    sum+=v
    data.append(v) # store the data
assert(math.isclose(sum,1))

h=0
for i in range(len(data)):
    h-= data[i]*math.log2(data[i])
print(data)
print(h)


"""
python3 entropy.py nucleotides.txt
1.846
"""
