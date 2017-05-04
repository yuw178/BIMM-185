#!/usr/bin/python
from Bio import SeqIO
import gzip
ncbids = {}
outdata = []
with gzip.open('uniprot_sprot_archaea.dat.gz', 'rb') as infile:
    with open('outfile2.txt','w') as outfile:
        outdata.append(['NCBI taxid','Organism','Taxonomy'])
        for record in SeqIO.parse(infile,"swiss"):
            if str(record.annotations['ncbi_taxid'])[2:-2] not in ncbids:
                ncbids[str(record.annotations['ncbi_taxid'])[2:-2]] = 1
                ncbi = str(record.annotations['ncbi_taxid'])[2:-2]
                org = record.annotations['organism']
                tax = ', '.join(record.annotations['taxonomy'])
                outdata.append(map(str,[ncbi,org,tax]))
            else:
                print str(record.annotations['ncbi_taxid'])[2:-2]
        for line in outdata:
            outfile.write('\t'.join(map(str,line)) + '\n')
