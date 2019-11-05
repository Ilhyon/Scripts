#!/usr/bin/env python
# -*- coding: utf-8 -*-:

def rankBwt(bw):
    ''' Retourne les rangs '''
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots:
            tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots


def firstCol(tots):
    ''' retourne la premiere colonne '''
    first = {}
    totc = 0
    for c, count in sorted(tots.iteritems()):
        first[c] = (totc, totc + count)
        totc += count
    return first

def reverseBwt(bw):
    ''' Retourne le texte original de la transformation bw '''
    ranks, tots = rankBwt(bw)
    first = firstCol(tots)
    rowi = 0
    t = "$"
    while bw[rowi] != '$':
        c = bw[rowi]
        t = c + t
        rowi = first[c][0] + ranks[rowi]
    return t

def main():
    seqBWT1 = 'CGGATTAC$ATG'
    seqBWT2 = 'CT$GAGGTCTAA'
    exp = 'ebn$naa'
    print(reverseBwt(seqBWT1))
    print(reverseBwt(seqBWT2))
    print(reverseBwt(exp))

main()
