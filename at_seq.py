#!/usr/bin/env python3

import random
#random.seed(1) # comment-out this line to change sequence each time

# Write a program that stores random DNA sequence in a string
# The sequence should be 30 nt long
# On average, the sequence should be 60% AT
# Calculate the actual AT fraction while generating the sequence
# Report the length, AT fraction, and sequence
length = 30
dna=""
count=0
alphabet=3*"A"+2*"C"+3*"T"+ 2*"G"

for i in range(length):
    nt=random.choice(alphabet)
    if nt=="A" or nt=="T":
        count+=1
        dna+=nt
print(length,count/length,dna)

"""
30 0.6666666666666666 ATTACCGTAATCTACTATTAAGTCACAACC
"""
