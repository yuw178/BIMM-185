#!/usr/bin/env python

with open('blast1_atvec.out') as infile: 
    alllines = [string.strip().split('\t') for string in infile.readlines()]
    for line in alllines:
        line = line.append(str(int(line[8])/int(line[3])))
with open('blast1_atvec.out','w') as outfile:
    for line in alllines:
        outfile.write('\t'.join(map(str,line))+'\n')
