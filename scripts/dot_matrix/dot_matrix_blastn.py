#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

usage: python %prog multi_fasta_1.fa multi_fasta_2.fa
version: Python 3.7
'''

import os,sys,glob,io,subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

matplotlib.rcParams['lines.linewidth']=1.5
matplotlib.rcParams['axes.linewidth']=1.5
matplotlib.rcParams['xtick.major.width']=1.5
matplotlib.rcParams['ytick.major.width']=1.5
matplotlib.rcParams['font.size']=15


# parameters for blast
evalue='1e-03'
word_size=7
thread_num=4
dust_yn='no'

# files
f_path_1=sys.argv[1]
f_path_2=sys.argv[2]
blastdb_name='tmp_blastdb_' + os.path.splitext(os.path.basename(f_path_2))[0]
title=os.path.basename(f_path_1) +' vs '+ os.path.basename(f_path_2)


# defs
def parse_fasta(path_to_file):
    tmp={}
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=seq
                header=line.strip().replace('>', '')
                seq=''
            elif '>' in line and not seq:
                header=line.strip().replace('>', '')
            else:
                seq += line.strip()
        tmp[header]=seq
    return tmp

def calc_fasta_len(fasta_d):
    tmp=[]
    for i in fasta_d:
        tmp.append(len(fasta_d[i]))
    return tmp

def parse_blast_result(dict):
    out=[]
    for i in dict:
        matched=[]  # match
        n=0
        for c in dict[i]['match']:
            if c == '|':
                matched.append(n)
            n += 1
        query_pos=[]  # query
        n=0
        m=0
        for c in dict[i]['query']:
            if m in matched: query_pos.append(n)
            if not c == '-': n += 1
            m += 1
        subj_pos=[]  # subj
        n=0
        m=0
        for c in dict[i]['sbjct']:
            if m in matched: subj_pos.append(n)
            if not c == '-': n += 1
            m += 1
        match_pos=[]
        start=''
        for n in range(len(matched)):  # reshape
            if not start and not start == 0:
                start=n
            elif not matched[n]+1 in matched:
                end=n
                match_pos.append([start, end])
                start=''
        for start,end in match_pos:
            if dict[i]['frame'] >= 1:
                out.append([dict[i]['qstart'] + query_pos[start],
                            dict[i]['qstart'] + query_pos[end],
                            dict[i]['sstart'] + subj_pos[start],
                            dict[i]['sstart'] + subj_pos[end]])
            else:
                out.append([dict[i]['qstart'] + query_pos[start],
                            dict[i]['qstart'] + query_pos[end],
                            dict[i]['sstart'] - subj_pos[start],
                            dict[i]['sstart'] - subj_pos[end]])
    return out


# parse fasta
fa1=parse_fasta(f_path_1)
fa2=parse_fasta(f_path_2)

# calc. fasta length
fa_len1=calc_fasta_len(fa1)
fa_len2=calc_fasta_len(fa2)

# number of fasta seqs
seq_num1=len(fa_len1)
seq_num2=len(fa_len2)


# make blastdb
cmd='makeblastdb -in %s -out %s -dbtype nucl -parse_seqids' % (f_path_2, blastdb_name)
process=subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
process.wait()


# BLASTn
res={}
out=NcbiblastnCommandline(query=f_path_1, db=blastdb_name, evalue=evalue, word_size=word_size,
                           num_threads=thread_num, outfmt=5, dust=dust_yn)()[0] # xml
blast_records=NCBIXML.parse(io.StringIO(out))
for blast_record in blast_records:
    query=blast_record.query.replace(' No definition line', '')
    res[query]={}
    for alignment in blast_record.alignments:
        subj=alignment.title.replace(' No definition line', '')
        tmp={}
        n=0
        for hsp in alignment.hsps:              # Bio.Blast.Record.HSP
            tmp[n]={}
            tmp[n]['qstart']=hsp.query_start - 1  # blast result = 1-start
            tmp[n]['qend']=  hsp.query_end
            tmp[n]['sstart']=hsp.sbjct_start
            tmp[n]['send']=  hsp.sbjct_end - 1  # blast result = 1-start
            tmp[n]['frame']= hsp.frame[1]
            tmp[n]['query']= hsp.query
            tmp[n]['match']= hsp.match
            tmp[n]['sbjct']= hsp.sbjct
            n += 1
        res[query][subj]=tmp


# prepare plt
y_len= 5 * seq_num1      # fig size
x_len= y_len * ( sum(fa_len1) / sum(fa_len2) )

plt.figure(figsize=(x_len, y_len))

if (seq_num1 >= 2) and (seq_num2 >= 2):
    gs=gridspec.GridSpec(seq_num2,seq_num1, width_ratios=fa_len1, height_ratios=fa_len2)
elif seq_num1 >= 2:
    gs=gridspec.GridSpec(seq_num2,seq_num1, width_ratios=fa_len1)
elif seq_num2 >= 2:
    gs=gridspec.GridSpec(seq_num2,seq_num1, height_ratios=fa_len2)
else:
    gs=gridspec.GridSpec(seq_num2,seq_num1)

gs.update(hspace=0.05, wspace=0.15, left=0.05, right=0.95, bottom=0.1, top=0.95) # left=0.1, right=0.95, bottom=0.1, top=0.95


# plot
r=0
for subj in fa2:
    c=0
    for query in fa1:
        ax=plt.subplot(gs[c+(r*seq_num1)])
        ax.set_xlim(0,fa_len1[c])
        ax.set_ylim(0,fa_len2[r])
        if c == 0:
            ax.set_ylabel(subj)
        if r == (seq_num2 - 1):
            ax.set_xlabel(query)
        if (c == 0) and (r < (seq_num2 - 1)):
            ax.xaxis.set_ticklabels([])
        elif (c > 0) and (r < (seq_num2 - 1)):
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
        elif (c > 0) and (r == (seq_num2 - 1)):
            ax.yaxis.set_ticklabels([])
        
        if query in res:
            if subj in res[query]:  # plot blast result
                df=pd.DataFrame(parse_blast_result(res[query][subj]))
                df=df.values
                ax.plot(df[:,0:2].T, df[:,2:4].T, color='black', alpha=0.5, linewidth=1.5) # dot plot line conf; , alpha=0.5
            else:  # plot 'no_hit'
                plt.text(fa_len1[c]/2, fa_len2[r]/2, 'No_hit', alpha=0.3, ha='center', va='center')
        else:  # plot 'no_hit'
            plt.text(fa_len1[c]/2, fa_len2[r]/2, 'No_hit', alpha=0.3, ha='center', va='center')
        c += 1
    r += 1

plt.suptitle(title)
base=os.path.splitext(os.path.basename(f_path_2))[0]
plt.savefig('plot_out_%s.pdf' % base)


# delete blastdb
blast_db_files=glob.glob('%s.n*' % blastdb_name)
if len(blast_db_files) >= 1:
    for i in blast_db_files:
        if os.path.isfile(i):
            os.remove(i)

