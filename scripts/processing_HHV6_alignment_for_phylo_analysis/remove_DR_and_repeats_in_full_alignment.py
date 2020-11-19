#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

usage: python %prog alignment.fa [HHV6A, HHV6B]
version: Python 3.7
'''


import os,sys,glob


f=sys.argv[1]
A_or_B=sys.argv[2]  # HHV6A or HHV6B


a_rep=[[0, 8089], [127548, 128233], [131076, 131854], [140075, 140951], [151288, 159378]]  # annotation in 'NC_001664.4'; DRL, R1, R2, R3, DRR
b_rep=[[0, 8793], [9314, 9510], [129045, 129681], [133500, 133863], [133981, 134076], [140081, 142691], [153321, 162114]]  # annotation in 'NC_000898.1'; DRL, R0, R1, R2A, R2B, R3, DRR


a_ref='>X83413.2 Human betaherpesvirus 6A, variant A DNA, complete virion genome, isolate U1102'   # specify the reference HHV-6A fasta header
b_ref='>AF157706.1 Human herpesvirus 6B strain Z29, complete genome'   # specify the reference HHV-6B fasta header


def parse_fasta(path_to_file):
    tmp={}
    seq=[]
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=''.join(seq)
                header=line.strip()
                seq=[]
            elif '>' in line and not seq:
                header=line.strip()
            else:
                seq.append(line.strip())
        tmp[header]=''.join(seq)
    return tmp


def remove_rep(f, hhv_rep, hhv_ref):
    starts=set()
    ends_minus_one=set()
    prev_e=-1
    for s,e in hhv_rep:
        if prev_e > 0:
            starts.add(prev_e)
            ends_minus_one.add(s - 1)
        prev_e=e
    
    fa=parse_fasta(f)
    retain_align_starts=[]
    retain_align_ends=[]
    next_pos_is_end=False
    align_n=0
    genome_n=0
    for c in fa[hhv_ref]:
        if genome_n in starts:
            retain_align_starts.append(align_n)
        elif next_pos_is_end is True:
            if c == '-':
                pass
            else:
                retain_align_ends.append(align_n)
                next_pos_is_end=False
        elif genome_n in ends_minus_one:
            next_pos_is_end=True
        if not c == '-':
            genome_n += 1
        align_n += 1
    retain_align_pos=[]
    for s,e  in zip(retain_align_starts, retain_align_ends):
        retain_align_pos.append([s, e])
    for h in fa:
        tmp=[]
        for s,e in retain_align_pos:
            tmp.append(fa[h][s:e])
        seq=''.join(tmp)
        tmp=[]
        for i in range(0, len(seq), 60):
            tmp.append(seq[i:i+60])
        fa[h]='\n'.join(tmp)
    filename,_=os.path.splitext(f)
    with open('%s_rmrep.fa' % filename, 'w') as outfile:
        for h in fa:
            outfile.write('%s\n%s\n' % (h, fa[h]))


if A_or_B == 'HHV6A':
    remove_rep(f, a_rep, a_ref)
elif A_or_B == 'HHV6B':
    remove_rep(f, b_rep, b_ref)
