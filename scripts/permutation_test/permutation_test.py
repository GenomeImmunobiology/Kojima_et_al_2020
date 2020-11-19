#!/usr/bin/env python

'''
Copyright (c) 2020 RIKEN
All Rights Reserved
See file LICENSE for details.

usage: python %prog
version: Python 3.7
'''


import os,sys,math,subprocess
from statistics import mean


original_poss=['/path_to/GWAS_peaks.txt']    # original position; bed file
target_poss=['/path_to/all_piG_concatenated.bed',
            '/path_to/piG.hybrid.merged.bed',
            '/path_to/piG.pachytene.merged.bed',
            '/path_topiG.prepachytene.merged.bed']     # target position; bed file


chrom_size='/path_to/hg38.fa.fai'  # specify the human genome fasta index
len_extend=0  # extend +/- from original element
n_permut=1000

len_extend=str(len_extend)


def count_intersect_num(original_pos, target_pos, chrom_size, len_extend):
    cmd='bedtools intersect -a %s -b %s -wa | sort | uniq | wc -l' % (original_pos, target_pos)
    count=subprocess.check_output(cmd, shell=True).decode().strip()
    return int(count)

def count_randomized_intersect_num(original_pos, target_pos, chrom_size, len_extend):
    cmd='bedtools shuffle -noOverlapping -i %s -g %s | bedtools intersect -a "stdin" -b %s -wa | sort | uniq | wc -l'  % (original_pos, chrom_size, target_pos)
    count=subprocess.check_output(cmd, shell=True).decode().strip()
    return int(count)


if os.path.exists('out_permutation_test.txt') is False:
    with open('out_permutation_test.txt', 'w') as outfile:  # randomized intersect num result
        outfile.write('target_genomic_pos\tsubjects\telongation\tobserved_num\trandomized_num\tfold_enrichment\tpval\n')


permut_num_out={}

with open('out_permutation_test.txt', 'a') as outfile:  # pval result
    for original_pos in original_poss:
        for target_pos in target_poss:
            orig_count=count_intersect_num(original_pos, target_pos, chrom_size, len_extend)   # count intersected element num
            random_l=[]
            for i in range(n_permut):  # randomize
                count_random=count_randomized_intersect_num(original_pos, target_pos, chrom_size, len_extend)
                random_l.append(count_random)
            permut_num_out[target_pos]=random_l
            mean_random=mean(random_l)
            if mean_random > 0:
                enrich= (orig_count / mean_random)  # enrich
            else:
                enrich='NA'
            count_exceed_num=0
            for r in random_l:
                if r >= orig_count:
                    count_exceed_num += 1
            if count_exceed_num > 0:
                pval=str(count_exceed_num / n_permut)  # pval
            else:
                pval_min= (1 / n_permut)
                pval='<' + str(pval_min)
            tmp=[original_pos, target_pos, '(+/-)'+ str(len_extend) +'-nt', str(orig_count), str(mean_random), str(enrich), str(pval)]
            outfile.write('\t'.join(tmp) +'\n')
            print('\t'.join(tmp))


out=[]
for t in permut_num_out:
    l=[ str(i) for i in permut_num_out[t] ]
    out.append(t +'\t'+ '\t'.join(l))
with open('out_permutation_test.tsv', 'w') as outfile:  # randomized intersect num result
    outfile.write('\n'.join(out) +'\n')

