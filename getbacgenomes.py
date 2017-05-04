#!/usr/bin/python
import subprocess

with open('README') as infile:
    b = [line.strip().split() for line in infile.readlines()]
    bacteria = []
    for line in b:
        if len(line) == 9:
            bacteria.append(line)
    ids = bacteria[3:6]
    print ids
    for bac in ids:
        tocall = "wget -P " + bac[0] + " ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/" + bac[0] + "_*"
        print tocall
        p = subprocess.Popen(tocall, shell=True)
