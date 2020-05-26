#!/usr/bin/env python3

# Write a program that compares two files of names to find:
#	Names unique to file 1
#	Names unique to file 2
#	Names shared in both files


import sys

file1 = sys.argv[1]
file2 = sys.argv[2]


contents1 = {}
contents2 = {}

def dictionary(filename):
    dictionary = {}
    with open(filename) as fp:
        for word in fp.readlines():
            word = word.rstrip()
            dictionary[word] = True
    return dictionary
 
d1 = dictionary(file1)
d2 = dictionary(file2)

list1 = []
list2 = []
shared = []

for name in d1:
    if name in d2: shared.append(name)
    else:          list1.append(name)

for name in d2:
    if name not in d1: list2.append(name)

print(f'Unique to file1: {list1}')
print(f'Unique to file2: {list2}')
print(f'Shared in both files: {shared}')

"""
python3 checklist.py --file1 --file2
"""