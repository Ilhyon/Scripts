#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

import os
import argparse
import itertools
import pandas as pd
from Bio import SeqIO
from BCBio import GFF
from pprint import pprint

def overlaps(interval1, interval2):
    """Compute the distance or overlap between two interval.

    Compute : min(ends) - max(starts). If > 0, return the number of bp of
    overlap, if 0, they are book-ended and if < 0 return the distance in
    bp between them.

    :param interval1: contain the start and the end of a OQs.
	:type interval1: list of int
    :param interval2: contain the start and the end of a gff feature.
	:type interval2: list of int

    :returns: min(ends) - max(starts)
    :rtype: int
    """
    return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])

def overlapSeqFeature(intervalSeq, intervalFeat):
    """Get the overlap.

    If there is an overlap (>0) between an OQs and a feature, return the feature.

    :param intervalSeq: contain the start and the end of a OQs.
	:type intervalSeq: list of int
    :param intervalFeat: contain the start and the end of a gff feature.
	:type intervalFeat: list of int

    :returns: intervalFeat,  contain the start and the end of a gff feature.
    :rtype: list
    """
    overlap = overlaps(intervalSeq, intervalFeat)
    if overlap > 0:
        return intervalFeat

def importGff(gff):
    """Import a gf file into a dataframe.

    The file is almost not parsed, we only remove empty columns and we add one
    with a list containing the start and the end of the feature.

    :param gff: name of the gff file.
	:type gff: string

    :returns: df, each row is a feature and columns are Chromosome, Source,
        Feature, Start, End, Strand, Attributes, Coords.
    :rtype: dataFrame
    """
    try:
        df = pd.read_csv(gff, sep='\t', skiprows=23)
    except:
        print("This file couldn't be converted in data frame : " + gff)
    else:
        df.columns = ['Chromosome', 'Source', 'Feature','Start','End',
                    'NI1', 'Strand', 'NI2', 'Attributes']
        del(df['NI1'])
        del(df['NI2'])
        df = df[ df.Chromosome != '###']
        df = df.reset_index(drop=True)
        df['Coords'] = [ [df.Start[x], df.End[x]] for x in range(0,len(df))]
        return df

def getCommonElem(fasta, df):
    """Get OQs that fit to a gff feature.

    For each OQs, we search if there is an overlap with gff feature that are on
    the same chromosome and on the same strand.

    :param fasta: name of the fasta file containing OQs.
	:type fasta: string
    :param df: contains all feature of a gff file.
	:type df:dataFrame

    :returns: dicoCommon, {OQs name : [feature coord]}, the list of feature
        coord are only overlapping feature with the OQs key.
    :rtype: dictionary
    """
    dicoCommon = {}
    fastaOrigin = SeqIO.parse(open(fasta),'fasta')
    for f in fastaOrigin:
        name, sequence = f.id, str(f.seq)
        chr = name.split(':')[0].split('r')[1]
        coordSeq = name.split(':')[1].split('(')[0].split('-')
        coordSeq = [ int(coordSeq[0]), int(coordSeq[1]) ]
        strand = name.split('(')[1].split(')')[0]
        dfTmp = df[ df.Strand == strand ]
        dfTmp = dfTmp[ dfTmp.Chromosome == chr ]
        listCommon = [overlapSeqFeature(coordSeq, c) for c in dfTmp.Coords]
        listCommon = sorted([x for x in listCommon if x is not None])
        dicoCommon[name] = list(listCommon for listCommon,_ in itertools.groupby(listCommon))
    return dicoCommon

def getCoordDico(df):
    """ Parse a gff dataframe into a dictionary.

    To ease the matching between OQs and feature we make a dictionary containing
    feature Coords as key and a list of corresponding feature as values.

    :param df: contains all feature of a gff file.
	:type df:dataFrame

    :returns: dico, {tuple(coords) : [feature]}.
    :rtype: dictionary
    """
    dico = {}
    dfTmp = df[['Coords', 'Feature']]
    dfTmp.Coords = [tuple(x) for x in dfTmp.Coords]
    for index, row in dfTmp.iterrows():
        if tuple(row.Coords) not in dico:
            dico[tuple(row.Coords)] = []
        dico[tuple(row.Coords)].append(row.Feature)
    # pprint(dfTmp.set_index('Coords').T.to_dict('list'))
    return dico

def main(fasta, gff):
    """ Print OQs corresponding to features with feature coords and names.

    First data are imported and parsed (gff as dataframe and dictionary) and
    then common element are retrieve between Oqs and gff feature.

    :param fasta: name of the fasta file containing OQs.
	:type fasta: string
    :param gff: name of the gff file.
	:type gff: string
    """
    df = importGff(gff)
    dicoCoord = getCoordDico(df)
    dicoCommon = getCommonElem(fasta, df)
    for coordSeq in dicoCommon:
        [print(coordSeq, x, dicoCoord[tuple(x)]) for x in dicoCommon[coordSeq] if tuple(x) in dicoCoord]

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'readBlast')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    gffFile = path+'Saccharomyces_cerevisiae.R64-1-1.87.gff3'
    # fastaFile = path+'GSM3003553_Saccaromyces_all_w15_th-1_minus.hits.max.K.w50.35.fasta'
    fastaFile = path+'GSM3003553_Saccaromyces_all_w15_th-1_plus.hits.max.K.w50.35.fasta'
    main(fastaFile, gffFile)
