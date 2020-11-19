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


a_rep=[[0, 376], [1682, 1730], [2302, 2367], [2369, 2451], [2692, 2733], [3149, 3181], [3433, 3502], [3626, 3670], [7483, 7519], [7655, 500000]]  # 8089
b_rep=[[0, 393], [1926, 2011], [2674, 2717], [3013, 3067], [3670, 3713], [3959, 3988], [8248, 500000]]  # 8793


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
    removed=[]
    for h in fa:
        tmp=[]
        for s,e in retain_align_pos:
            tmp.append(fa[h][s:e])
        seq=''.join(tmp)
        tmp=[]
        for i in range(0, len(seq), 60):
            tmp.append(seq[i:i+60])
        count_n=0
        total=0
        for partial in tmp:
            count_n += partial.count('N')
            count_n += partial.count('n')
            total += len(partial.replace('-', ''))
        if (count_n / total) > 1:
            removed.append(h)
        else:
            fa[h]='\n'.join(tmp)
    filename,_=os.path.splitext(f)
    with open('%s_rmrep.fa' % filename, 'w') as outfile:
        for h in fa:
            outfile.write('%s\n%s\n' % (h, fa[h]))
    print(';'.join(removed))
    print('removed due to high N = %d' % len(removed))


if A_or_B == 'HHV6A':
    remove_rep(f, a_rep, a_ref)
elif A_or_B == 'HHV6B':
    remove_rep(f, b_rep, b_ref)
