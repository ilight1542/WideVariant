#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 12:13:33 2022

@author: evanqu
"""
import numpy as np
import os
import gzip
import sys
import argparse
from Bio.Data import IUPACData
import gus_helper_functions as ghf
from math import floor
from gus_helper_functions import round_half_up, convert_chrpos_to_abspos


def get_fq_score_of_simple_call(vcf_ref, vcf_alt, vcf_info, remove_ambigous_call = True):
    ## only consider simple calls from vcf (alt is present; no multiple calls in alt, ref&alt same length; ref=1nt)
    unambiguous_nts = ['A', 'T', 'C', 'G']
    alt_not_multiallelic = ("," not in vcf_alt)
    alt_ref_same_length = (len(vcf_alt) == len(vcf_ref))
    ref_1nt = (len(vcf_ref)==1)
    if remove_ambigous_call:
        ambiguous_nts = [nt for nt in IUPACData.ambiguous_dna_values.keys() if nt not in unambiguous_nts]
        ref_unambiguous = all([ref not in ambiguous_nts for ref in vcf_ref])
        alt_unambiguous = all([alt not in ambiguous_nts for alt in vcf_alt])
        keep_pos = (ref_unambiguous & alt_unambiguous)
    else:
        keep_pos = True
    if (vcf_alt) and alt_not_multiallelic and alt_ref_same_length and ref_1nt and keep_pos:
        #find and parse quality score
        vcf_info_l = vcf_info.split(';')
        entrywithFQ_idx=[x for x in vcf_info_l if x.startswith('FQ')]
        if (entrywithFQ_idx == []) or len(entrywithFQ_idx) > 1: ## non ore more than one entry identified
            print('No or more than 1 FQ entries are found in provided vcf')
            print(f'The vcf info field was: {vcf_info}')
            sys.exit(0)
        else:
            entrywithFQ_str = entrywithFQ_idx[0]
            fq=entrywithFQ_str[entrywithFQ_str.index("=")+1:]
            return float(fq)
    else:
        return np.nan

def vcf_to_quals_snakemake(path_to_vcf_file,output_path_to_quals,REFGENOMEDIRECTORY):
    '''Python version of vcf_to_quals_snakemake.py
    Given a vcf file with one file per line, grabs FQ score for each positions. Ignores lines corresponding to indels

    Args:
        path_to_vcf_file (str): Path to .vcf file.
        output_path_to_quals (str): Path to output quals file
        REFGENOMEDIRECTORY (str): Path to reference genome directory.

    Returns:
        None.
        
    '''
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(REFGENOMEDIRECTORY)
    
    #initialize vector to record quals
    quals = np.zeros((genome_length,1), dtype=int)
    
    print(f"Loaded: {path_to_vcf_file}")
    
    with gzip.open(path_to_vcf_file, 'rt') as fid:
        for line in fid:
            if not line.startswith("#"):
                lineinfo = line.strip().split('\t')
                
                ## convert to position in genome from chromosome and position on chromosome
                chromo = lineinfo[0]
                chr_idx = np.where(chromo == scaf_names)[0][0] +1 ## +1 as chrpos2index uses 1-indexed chromosomes
                position_on_chr = int(lineinfo[1]) #1-indexed
                chr_pos_array = np.array( ([chr_idx, position_on_chr], ) ) ## 2D-np array required!
                position = ghf.chrpos2index(chr_pos_array, chr_starts)
                #position = ghf.convert_chrpos_to_abspos(chromo, position_on_chr, chr_starts, scaf_names)

                ## extract quality score if an alternative allele called
                ref_nt=lineinfo[3]
                alt_nt=lineinfo[4]
                
                #only consider for simple calls (not indel, not ambiguous)
                info_col = lineinfo[7]
                #only consider for simple calls (not indel, not ambiguous)
                fq_score = get_fq_score_of_simple_call(ref_nt, alt_nt, info_col)
                #If already a position wiht a stronger FQ here, don;t include this
                #More negative is stronger
                if fq_score < quals[position-1]:
                    quals[position-1]= int(round_half_up(fq_score)) #python int(fq) will by default round down, round matches matlab behavior; -1 important to convert position (1-indexed) to python index
        
    #save
    outpath = os.getcwd() + '/' + output_path_to_quals
    np.savez_compressed(outpath, quals=quals)
        
    print(f"Saved: {output_path_to_quals}")
    
    return

#%%
if __name__ == '__main__':
    
    SCRIPTS_DIR="scripts"
    sys.path.insert(0, SCRIPTS_DIR)
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', type=str, help='Path to input vcf file',required=True)
    parser.add_argument('-r', type=str, help='Path to reference genome directory',required=True)
    parser.add_argument('-o', type=str, help='Path to output quals file (.npz)', required=True)
    
    args = parser.parse_args()
    
    vcf_to_quals_snakemake(args.i,args.o,args.r)
