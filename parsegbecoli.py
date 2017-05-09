#!/usr/bin/env python

from Bio import SeqIO
import gzip
outdata = []
with gzip.open('GCF_000005845.2_ASM584v2_genomic.gbff.gz', 'rb') as infile:
    with open('../../genes.txt','w') as outfile:
        #indexing
        i = 1
        j = 1
        k = 1
        l = 1
        genes = []
        exons = []
        synonyms = []
        references = []
        functions = []
        #using seqIO to parse
        for record in SeqIO.parse(infile,"genbank"):
            #getting genome info
            taxid = int(str(record.features[0].qualifiers['db_xref'])[8:-2])
            short = record.name
            longnm = record.description
            sz = len(record.seq)
            domain = record.annotations['taxonomy'][0]
            accession = record.name
            date = record.annotations['date']
            genomeline = [i, taxid, short, longnm, sz, domain, accession, date]
            #counting genes
            cds = 0
            for feature in record.features:
                if feature.type == "CDS":
                    cds += 1
                    #getting gene info
                    if 'locus_tag' in feature.qualifiers:
                        locus = str(feature.qualifiers['locus_tag'])[2:-2]
                    else:
                        locus = '-'
                    if 'gene' in feature.qualifiers:
                        name = str(feature.qualifiers['gene'])[2:-2]
                    else:
                        name = '-'
                    strand = str(feature.strand)
                    #getting start/stop/length info
                    numexons = len(str(feature.location).split(','))
                    if numexons > 1:
                        coords = str(feature.location)[5:-1].split(', ')
                    else:
                        coords = [str(feature.location)]
                    bp = 0
                    #getting exon info
                    for item in coords:
                        start = int(item[1:-4].split(':')[0])
                        end = int(item[1:-4].split(':')[1])
                        length = end-start
                        bp += length
                        exonline = [k,l,start,end,length]
                        exons.append(exonline)
                        l+=1
                        #more gene info
                    if 'product' in feature.qualifiers:
                        protein = str(feature.qualifiers['product'])[2:-2]
                    else:
                        protein = '-'
                    geneline = [k,i,j,locus,name,strand,numexons,bp,protein]
                    genes.append(geneline)
                    #getting synonyms, references, functions
                    if 'gene_synonym' in feature.qualifiers:
                        syns = str(feature.qualifiers['gene_synonym'])[2:-2]
                        synonyms.append([k,syns])

                    if 'db_xref' in feature.qualifiers:
                        refs = str(feature.qualifiers['db_xref'])[1:-1].split(', ')
                        if 'protein_id' in feature.qualifiers:
                            references.append([k,'refseq', str(feature.qualifiers['protein_id'])[2:-2]])
                        for ref in refs:
                            db = ref[1:-1].split(':')[0]
                            dbid = ref[1:-1].split(':')[1]
                            references.append([k,db,dbid])

                    if 'function' in feature.qualifiers:
                        function = str(feature.qualifiers['function'])[2:-2]
                        functions.append([k,function])
                    k+=1
            structure = record.annotations['topology']
            if structure in ['circular', 'linear']:
                reptype = 'chromosome'
            else:
                reptype = 'plasmid'
            repliconline = [j,i,record.description,cds,len(record.seq),reptype,structure,record.name,record.annotations['date']]
            print '\t'.join(map(str,repliconline))
            #writing to files
            for gene in genes:
                outfile.write('\t'.join(map(str,gene)) + '\n')
            with open('../../exons.txt','w') as outfile2:
                for exon in exons:
                    outfile2.write('\t'.join(map(str,exon)) + '\n')
            with open('../../synonyms.txt','w') as outfile3:
                for synonym in synonyms:
                    outfile3.write('\t'.join(map(str,synonym)) + '\n')
            with open('../../references.txt','w') as outfile4:
                for reference in references:
                    outfile4.write('\t'.join(map(str,reference)) + '\n')
            with open('../../functions.txt','w') as outfile5:
                for function in functions:
                    outfile5.write('\t'.join(map(str,function)) + '\n')
