#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import sys
from Bio import SeqIO

def split_fasta (fastafile, dirname):

    infile = open(fastafile)
    for record in SeqIO.parse(infile, "fasta"):
		idSeq = record.id
		seq = str(record.seq)
		filename = record.id.split("|")[0]
		outfile = open ("%s/%s.fasta" % (dirname, filename), "w")
		outfile.write(">"+idSeq+"\n")
		outfile.write(seq)
		outfile.close()    
        

split_fasta(sys.argv[1], sys.argv[2])
