#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:36:31 2022

@author: evanqu
"""
import numpy as np
import os
import gzip
import sys
import argparse
from Bio.Data import IUPACData
import gus_helper_functions as ghf

#%%
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

def generate_positions_single_sample(path_to_variant_vcf,path_to_output_positions,maxFQ,REFGENOMEDIRECTORY,outgroup_bool):
    '''Python version of generate_positions_single_sample_snakemake.m

    Args:
        path_to_variant_vcf (str): Path to .variant.vcf.gz file.
        path_to_output_positions (str): Output path to positions file (.pickle.gz)
        maxFQ (int): Purity threshold for including position.
        REFGENOMEDIRECTORY (str): Path to reference genome directory.
        outgroup_bool (bool): Whether this sample is outgroup or not.

    Returns:
        None.

    '''    
    print(f"Currently examining the following vcf file: {path_to_variant_vcf}\n")
    print(f"FQ threshold: {int(maxFQ)}")
    
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(REFGENOMEDIRECTORY)

    # Initialize boolean vector for positions to include as candidate SNPs that
    # vary from the reference genome
    include = np.zeros((genome_length,1))
    
    #For outgroup samples only
    if outgroup_bool:
        Var_positions = ghf.p2chrpos(np.nonzero(include)[0]+1,chr_starts)
        np.savez_compressed(path_to_output_positions, Positions = Var_positions)
        print("Outgroup sample - no positions collected")
        return
    
    with gzip.open(path_to_variant_vcf, 'rt') as fid:
        for line in fid:
            if not line.startswith("#"):
                lineinfo = line.strip().split('\t')
                
                chromo = lineinfo[0]
                chr_idx = np.where(chromo == scaf_names)[0][0] +1 ## +1 as chrpos2index uses 1-indexed chromosomes
                position_on_chr = int(lineinfo[1]) #1-indexed
                chr_pos_array = np.array( ([chr_idx, position_on_chr], ) ) ## 2D-np array required!
                position = ghf.chrpos2index(chr_pos_array, chr_starts)
                #position = ghf.convert_chrpos_to_abspos(chromo, position_on_chr, chr_starts, scaf_names)
                    
                alt_nt = lineinfo[4]
                ref_nt = lineinfo[3]
                info_col = lineinfo[7]
                #only consider for simple calls (not indel, not ambiguous)
                fq_score = get_fq_score_of_simple_call(ref_nt, alt_nt, info_col)
                    
                if fq_score < maxFQ: #better than maxFQ
                    include[position-1] = 1 #-1 converts position (1-indexed) to index
    
    #+1 converts index back to position for p2chrpos
    Var_positions = ghf.p2chrpos(np.nonzero(include)[0]+1,chr_starts)
    
    #save
    outpath = os.getcwd() + '/' + path_to_output_positions
    np.savez_compressed(outpath, Positions = Var_positions)
        
    print(f"{len(Var_positions)} variable positions found passing quality threshold")
    
    return


if __name__ == '__main__':
    
    SCRIPTS_DIR="scripts"
    sys.path.insert(0, SCRIPTS_DIR)
    
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', type=str, help='Path to input variant vcf',required=True)
    parser.add_argument('-r', type=str, help='Path to reference genome directory',required=True)
    parser.add_argument('-o', type=str, help='Path to output positions file', required=True)
    parser.add_argument('-b', type=int, help='Outgroup boolean', required=True)
    parser.add_argument('-q', type=int, help='MaxFQ threshold', required=True)
    
    args = parser.parse_args()
    
    generate_positions_single_sample(args.i,args.o,args.q,args.r,args.b)

    
