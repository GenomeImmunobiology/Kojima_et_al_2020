#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

usage: python %prog dir_to_files
version: Python 3.7
'''

import os,sys,glob
from Bio.Blast import NCBIXML
import pysam
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

matplotlib.rcParams['lines.linewidth']=0.5
matplotlib.rcParams['axes.linewidth']=0.5
matplotlib.rcParams['xtick.major.width']=0.5
matplotlib.rcParams['ytick.major.width']=0.5
matplotlib.rcParams['font.size']=5


# virus fai, smrv
virus=[0, 8785]
lltr=[0, 456]
rltr=[8329, 8785]

near_len=10000

# load hu genome fai
fai={}
f='/path_to/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai'      # specify path to the reference human genome index
with open(f) as infile:
    for line in infile:
        ls=line.split()
        fai[ls[0]]= int(ls[1])
        if 'chrY' in line:
            break
chrs=list(fai.keys())
chrs_set=set(chrs)

chr_boxs={}
sum_len=0
spacer=10000000  # 10MB
for chr in fai:
    chr_boxs[chr]=[sum_len, fai[chr]]
    sum_len= sum_len + fai[chr] + spacer

chr_extend={}
for chr in chr_boxs:
    chr_extend[chr]=chr_boxs[chr][0]

out=[]

def plot(dir):
    try:
        # count all reads
        f='%s/mapped_to_target.bam' % dir   # specify a bam file with reads mapping to viruses
        mapped_rn=set()
        with pysam.AlignmentFile(f, 'rb') as infile:
            for line in infile:
                if line.is_secondary is False:
                    mapped_rn.add(line.query_name)
        
        # count hybrid reads
        f='%s/mapped_to_target_f8.bam' % dir   # specify a bam file with reads mapping to viruses, must only contain reads with sam flag 8
        hybrid_n=0
        with pysam.AlignmentFile(f, 'rb') as infile:
            for line in infile:
                if line.is_secondary is False:
                    hybrid_n += 1
        
        sample_id=dir.split('/')[2].split('_')[0]
        # load blastn results
        f='%s/blastn_mapped_to_target_f4_s.xml' % dir   # specify a blastn output file which contains results of read-mapping to the reference human genome, must be xml format.
        mapped={}
        blast_records=NCBIXML.parse(open(f))
        for blast_record in blast_records:
            if len(blast_record.alignments) == 1:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:              # Bio.Blast.Record.HSP
                        chr=alignment.title.split()[0]
                        if chr in chrs_set:
                            read_name=blast_record.query.split('/')[0]
                            chr_pos=hsp.sbjct_end + chr_extend[chr]
                            if hsp.strand[1] == 'Plus':
                                chr_strand='+'
                            else:
                                chr_strand='-'
                            mapped[read_name]=[chr_pos, chr_strand]

        # load mates; virus side
        f='%s/mapped_to_target_f8.bam' % dir   # specify a bam file with reads mapping to viruses, must only contain reads with sam flag 8
        with pysam.AlignmentFile(f, 'rb') as infile:
            for line in infile:
                if line.query_name in mapped:
                    if line.is_secondary is False:
                        if line.is_reverse is False:
                            virus_strand='+'
                            pos=line.reference_end
                        else:
                            virus_strand='-'
                            pos=line.reference_start
                        if rltr[0] <= pos <= rltr[1]:
                            pos= pos - rltr[0]
                        mapped[line.query_name].append(pos)
                        mapped[line.query_name].append(virus_strand)

        for read in mapped:
            if not len(mapped[read]) == 4:
                print(read)
        
        # identify near breakpoints
        breakpoints=set()
        for read in mapped:
            chr_pos,chr_strand,v_pos,v_strand=mapped[read]
            for r in mapped:
                if not r == read:
                    chr_pos2,chr_strand2,v_pos2,v_strand2=mapped[r]
                    if -1 * near_len < (chr_pos - chr_pos2) < near_len:
                        if not chr_strand == chr_strand2:
                            breakpoints.add(min([chr_pos, chr_pos2]))
        info='%s\t%d\t%d\t%d\t%d\n' % (sample_id, len(mapped_rn), hybrid_n, len(mapped), len(breakpoints))
        out.append(info)
        print(info, end='')
        
        # plot
        ymax=2
        box_height=0.1

        fig=plt.figure(figsize=(4, 2))  # (x, y)
        ax=fig.add_subplot(111)

        # draw chr
        a=0.25
        for chr in chr_boxs:
            a=0.25 if a == 0.5 else 0.5
            rect=matplotlib.patches.Rectangle((chr_boxs[chr][0], 0), chr_boxs[chr][1], box_height, color='k', alpha=a, ec=None)
            ax.add_patch(rect)
        xmax= chr_boxs[chr][0] + chr_boxs[chr][1]
        vpos_coeff= np.floor(xmax / virus[1])

        # draw virus
        rect=matplotlib.patches.Rectangle((virus[0] * vpos_coeff, ymax - box_height), (virus[1] - virus[0]) * vpos_coeff, box_height, color='tab:gray', alpha=0.25, ec=None)
        ax.add_patch(rect)
        rect=matplotlib.patches.Rectangle((lltr[0] * vpos_coeff, ymax - box_height), (lltr[1] - lltr[0]) * vpos_coeff, box_height, color='tab:gray', alpha=0.5, ec=None)
        ax.add_patch(rect)
        rect=matplotlib.patches.Rectangle((rltr[0] * vpos_coeff, ymax - box_height), (rltr[1] - rltr[0]) * vpos_coeff, box_height, color='tab:gray', alpha=0.5, ec=None)
        ax.add_patch(rect)

        scatter_pos_x,scatter_pos_y,scatter_neg_x,scatter_neg_y=[],[],[],[]
        for read in mapped:
            chr_pos,chr_strand,v_pos,v_strand=mapped[read]
            if chr_strand == '+':
                scatter_pos_x.append(chr_pos)
                scatter_pos_y.append(box_height)
                y_pos1=box_height
            else:
                scatter_neg_x.append(chr_pos)
                scatter_neg_y.append(0)
                y_pos1=0
            if v_strand == '+':
                scatter_pos_x.append(v_pos * vpos_coeff)
                scatter_pos_y.append(ymax)
                y_pos2=ymax
            else:
                scatter_neg_x.append(v_pos * vpos_coeff)
                scatter_neg_y.append(ymax - box_height)
                y_pos2=ymax - box_height
            ax.plot([chr_pos, v_pos * vpos_coeff], [y_pos1, y_pos2], color='tab:gray', linewidth=0.25, alpha=0.25)

        ax.scatter(scatter_pos_x, scatter_pos_y, s=5, c='dodgerblue', alpha=0.25)
        ax.scatter(scatter_neg_x, scatter_neg_y, s=5, c='orangered', alpha=0.25)
        
        for pos in breakpoints:
            ax.plot([pos, pos], [0, box_height], color='k', linewidth=0.5, alpha=0.75)
        
        plt.suptitle(sample_id)
        plt.savefig('%s/plot_%s_smrv_hybrid.pdf' % (dir, sample_id))
        plt.close()
    except:
        info='%s\t%d\t%d\tNA\tNA\n' % (sample_id, len(mapped_rn), hybrid_n)
        out.append(info)
        print(info, end='')


dir=sys.argv[1]
plot(dir)

with open('1kGP_smrv_hybrid_summary.txt', 'w') as outfile:
    outfile.write(''.join(out))
