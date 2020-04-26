#!/usr/bin/env python3

from math import sqrt
import fileinput

# Write a program that computes typical stats
# Count, Min, Max, Mean, Std. Dev, Median
# No, you cannot import any other modules!
count = 0
min = 0
max = 0
mean = 0
std_sum = 0
med = 0
sum = 0
sd = 0
data = []

#Loop that finds Counts and Sum
for line in fileinput.input():
    if line.startswith("#"): continue
    value = line.rstrip()
    value = float(value)
    sum += value
    count += 1
    data.append(value)

data.sort()
mean = sum/count
min  = data[0]
max  = data[-1]

if count % 2 == 0:
    med = (data[int(count/2 - 1)] + data[int(count/2)] ) /2
else:
    med = data[int((count+1)/2 -1)]

for i in range(len(data)):
    std_sum += ((data[i]-mean)**2)
sd = sqrt(std_sum/(count))


print(f'Count: {count}')
print(f'Minimum: {min}')
print(f'Maximum: {max}')
print(f'Mean: {mean}')
print(f'Std. dev: {sd:.3f}')
print(f'Median: {med}')

"""
python3 stats.py numbers.txt
Count: 10
Minimum: -1.0
Maximum: 256.0
Mean: 29.147789999999997
Std. dev: 75.777
Median 2.35914
"""
