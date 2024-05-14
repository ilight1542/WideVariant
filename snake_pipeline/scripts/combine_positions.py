#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 20:00:41 2022

@author: evanqu
"""
import numpy as np
import pickle
import os
import argparse
import gzip
import gus_helper_functions as ghf

#%%
def chrpos2index(chrpos,chr_starts):
    '''Python version of chrpos2index.m

    Args:
        chrpos (arr): px2 array of position and chromsome idx.
        chr_starts (arr): Vector of chromosome starts (begins at 0).

    Returns:
        p (arr): Vector of position indexes.

    '''
    if np.size(chrpos,0) < np.size(chrpos,1):
        chrpos=chrpos.T
        print('Reversed orientation of chrpos')
        
    if len(chr_starts) == 1:
        p=chrpos[:,1]
    else:
        p=chr_starts[chrpos[:,0]-1]+chrpos[:,1]

    return p

def generate_positions_snakemake(positions_files_list, REFGENOMEDIRECTORY):
    '''Python version of generate_positions_snakemake.m
    
    Args:
        paths_to_input_p_files (list): List of input positions files.
        REFGENOMEDIRECTORY (str): Path to reference genome.

    Returns:
        combined_pos (arr): Vector of variable positions across samples.

    '''
    
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(REFGENOMEDIRECTORY)
    
    # initialize vector to count occurrances of variants across samples
    timesvariant = np.zeros((genome_length,1))
    
    for path in positions_files_list:
        #load in positions array for sample
        tmp=np.load(os.getcwd() + '/' + path.rstrip('\n'))
        positions=tmp['Positions']
        
        if len(positions)>=1:
            x=chrpos2index(positions,chr_starts)
            
            timesvariant[x]=timesvariant[x]+1
    
    
    #Keep positions that vary from the reference in at least one sample but
    #that don't vary from the reference in ALL samples
    combined_pos = np.where((timesvariant > 0) & (timesvariant < len(positions_files_list)))[0]
    
    return combined_pos
    

def combine_positions(path_to_positions_files, path_to_output_p_file, path_to_outgroup_boolean_file, REFGENOMEDIRECTORY):
    
    #in_outgroup: booleans giving whether or not each sample is in the outgroup
    print("Processing outgroup booleans...")
    in_outgroup=[]
    with open(path_to_outgroup_boolean_file) as file:
        for line in file:
            in_outgroup.append(int(line.strip()))
    
    #Get positions on reference genome
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(REFGENOMEDIRECTORY)
    
    #Find positions with at least 1 fixed mutation relative to reference genome
    print('\n\nFinding positions with at least 1 fixed mutation...\n')
    
    positions_files_ls=[]
    with open(path_to_positions_files) as file:
        for line in file:
            positions_files_ls.append(line)
    ## include only ingroup samples (in_outgroup == 0)
    ingroup_files_to_include = np.array(in_outgroup) == 0
    positions_files_ingroup_ls = np.array(positions_files_ls)[ingroup_files_to_include]

    # print(f"\nIngroup paths used to generate positions: {positions_files_ls[include]}")
    cp = generate_positions_snakemake(positions_files_ingroup_ls,REFGENOMEDIRECTORY)
    print(f"Found {len(cp)} positions where provided vcfs called a fixed variant in at least one in-group sample \n")

    #Todo: Add candidate positions manually
    #Combine different types of positions
    allp=cp
    p=np.sort(allp)

    print("Saving list of all positions...")
    np.savez_compressed(os.getcwd() + '/' + path_to_output_p_file, p=p)
    
    return

#%%
if __name__ == '__main__':
    
    # SCRIPTS_DIR="scripts"
    # sys.path.insert(0, SCRIPTS_DIR)
        
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', type=str, help='Path to List of input positions files',required=True)
    parser.add_argument('-r', type=str, help='Path to reference genome directory',required=True)
    parser.add_argument('-o', type=str, help='Path to output positions file', required=True)
    parser.add_argument('-b', type=str, help='Outgroup boolean', required=True)
    
    args = parser.parse_args()
    
    print()

    combine_positions(args.i,args.o,args.b,args.r)
