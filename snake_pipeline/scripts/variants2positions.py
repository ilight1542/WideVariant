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
import gus_helper_functions as ghf

#%%
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
        Var_positions=ghf.p2chrpos(np.nonzero(include)[0],chr_starts)
        np.savez_compressed(path_to_output_positions, Positions = Var_positions)
        print("Outgroup sample - no positions collected")
        return
    
    with gzip.open(path_to_variant_vcf, 'rt') as fid:
        for line in fid:
            if not line.startswith("#"):
                lineinfo = line.strip().split('\t')
                
                chromo=lineinfo[0]
                position_on_chr=lineinfo[1] #1-indexed
                
                if len(chr_starts) == 1:
                    position=int(lineinfo[1])
                else:
                    if chromo not in scaf_names:
                        raise ValueError("Scaffold name in vcf file not found in reference")
                    position=int(chr_starts[np.where(chromo==scaf_names)]) + int(position_on_chr)
                    #chr_starts begins at 0
                    
                alt=lineinfo[4]
                ref=lineinfo[3]
                
                #only consider for simple calls (not indel, not ambiguous)
                if (alt) and ("," not in alt) and (len(alt) == len(ref)) and (len(ref)==1):
                    #find and parse quality score
                    xt_info = lineinfo[7].split(';')
                    entrywithFQ=[x for x in xt_info if x.startswith('FQ')]
                    entrywithFQ_str = entrywithFQ[0]
                    fq=entrywithFQ_str[entrywithFQ_str.index("=")+1:]
                    
                    if float(fq) < maxFQ: #better than maxFQ
                        include[position-1]=1
                        #-1 converts position (1-indexed) to index
    
    Var_positions=ghf.p2chrpos(np.nonzero(include)[0],chr_starts)
    
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

    
