#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
from Bio import AlignIO

filename = '/home/anais/Documents/Data/Alignment/Output/geneList_1000016'
# multiple_alignment = AlignIO.parse(filename, 'clustal')
alignment = AlignIO.read(open(filename), 'clustal')
print("Alignment length %i" % alignment.get_alignment_length())
# print(alignment)
for aliPos in range(0,alignment.get_alignment_length()):
    print(alignment[:, aliPos])
    print('-------')
# for record in alignment :
#     print(record.seq + " " + record.id)
#     print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
#             (record.annotations["start"],
#             record.strand,
#             record.srcSize,
#             record.size))
# for seqrec in multiple_alignment:
#         print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
#               (seqrec.annotations["start"],
#                seqrec.annotations["strand"],
#                seqrec.annotations["srcSize"],
#                seqrec.annotations["size"]))
