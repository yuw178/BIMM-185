#!/usr/bin/env python

with open('blast2_atvec.out') as infile:
    besthitsat = {}
    maxscore = {}
    scoresat = [line.strip().split('\t') for line in infile.readlines()]
    for line in scoresat:
        if line[0] not in besthitsat:
            besthitsat[line[0]] = line[1]
            maxscore[line[0]] = line[4]
        else:
            if float(maxscore[line[0]]) < float(line[4]):
                besthitsat[line[0]] = line[1]
                maxscore[line[0]] = line[4]

with open('blast2_ecvat.out') as infile:
    besthitsec = {}
    maxscore = {}
    scoresat = [line.strip().split('\t') for line in infile.readlines()]
    for line in scoresat:
        if line[0] not in besthitsat:
            besthitsec[line[0]] = line[1]
            maxscore[line[0]] = line[4]
        else:
            if float(maxscore[line[0]]) < float(line[4]):
                besthitsec[line[0]] = line[1]
                maxscore[line[0]] = line[4]
bidirectionalat = []
bidirectionalec = []
for item in besthitsat.keys():
    if besthitsat[item] in besthitsec:
        if besthitsec[besthitsat[item]] == item:
            bidirectionalat.append([item,besthitsat[item]])
            bidirectionalec.append([besthitsat[item],item])
with open('homology_at.txt','w') as outfile:
    for item in bidirectionalat:
        outfile.write('\t'.join(map(str,item)) + '\n')
with open('homology_ec.txt','w') as outfile:
    for item in bidirectionalec:
        outfile.write('\t'.join(map(str,item)) + '\n')
