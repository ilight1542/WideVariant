#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 22:14:06 2022

@author: evanqu
"""

import numpy as np
import sys
import gzip
import argparse
import gus_helper_functions as ghf
import pickle
from scipy.stats import ttest_ind, fisher_exact,ttest_1samp
from math import log10, floor

#%% Version history
#2022.02.08: Evan: Direct translation from pileup_to_diversity_matrix_snakemake.m
#2022.10.18, Arolyn: Now works when reference genome has lowercase letters or ambiguous letters
#2022.10.23, Arolyn: Updated comments on 40 statistics to have python indexing (0-39) as opposed to matlab indexing (1-40)

#%%Some notes

# This function saves the following 40 statistics for each position on the
# reference genome for the sample being analyzed:
# [A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd E I D]
# List of statistics by index:
# [0-3] A is the number of forward reads supporting A
# [4-7] a is the number of reverse reads supporting A
# [8-15] Aq is the average phred qualities of all A's
# [16-23] Am is the average mapping qualities of all A's
# [24-31] At is the average tail distance of all A's
# [32] -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# Ps is the p value for strand bias (fishers test)
# [33] -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# Pb is the p value for the base qualities being the same for the two
# different types of calls (1st major, 2nd major nt, either strand) (ttest)
# [34] -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# Pm is the p value for the mapping qualities being the same for the two
# different types of calls (ttest)
# [35] -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# Pftd is the p value for the tail distantces on the forward strand
# being the same for the two different types of calls (ttest)
# [36] -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# Pftd is the p value for the tail distantces on the reverse strand
# being the same for the two  different types of calls (ttest)
# [37] E is number of calls at ends of a read -- STATISTIC NOT COMPUTED OR RECORDED IN THIS VERSION -- 
# [38] I is number of reads supporting insertions in the +/- (indelregion) bp region
# [39] D is number of reads supporting deletions in the +/- (indelregion) bp region

# ChrStarts: is an array holding the indices in the position dimension
# corresponding to the start of a new chromsome.

def pileup2diversity(input_pileup, path_to_ref,min_reads_on_strand=20):
    """Grabs relevant allele info from mpileupfile and stores as a nice array 

    Args:
        input_pileup (str): Path to input pileup file.
        path_to_ref (str): Path to reference genome file
        
    """
    #Set parameters
    Phred_offset=33 #mpileup is told that the reads are in fastq format and 
                    #corrects so its always at Phred+33 when the mpileup comes out
    nts='ATCGatcg'
    nts_dict={nt: i for i, nt in enumerate(nts)}
    nts_ascii=[ord(nt) for nt in nts]
    num_fields=40
    indelregion=3 #region surrounding each p where indels recorded 
    #get reference genome + position information
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(path_to_ref)
    
    #init
    data = np.zeros((genome_length,num_fields),dtype=int) #format [[A T C G  a t c g],[...]]
    
    #####
    loading_bar=0

    #Read in mpileup file
    print(f"Reading input file: {input_pileup}")
    with open(input_pileup) as mpileup:
        for line in mpileup:
            loading_bar+=1
            if loading_bar % 50000 == 0:
                print('.', end = '')

            lineinfo = line.strip().split('\t')
            
            #holds info for each position before storing in data
            temp = np.zeros((num_fields),dtype=int)
            
            chromo = lineinfo[0]
            #position (absolute)
            if len(chr_starts) == 1:
                position=int(lineinfo[1])
            else:
                if chromo not in scaf_names:
                    raise ValueError("Scaffold name in pileup file not found in reference")
                position=int(chr_starts[np.where(chromo==scaf_names)]) + int(lineinfo[1])
                #chr_starts starts at 0
            
            #ref allele
            ref_str = lineinfo[2]; # reference allele from pileup (usually A T C or G but sometimes a different symbol if nucleotide is ambiguous)
            if ref_str in nts_dict.keys():
                ref=nts_dict[ref_str] # convert to 0123
                if ref >= 4:
                    ref = ref - 4
            else:
                ref=None # for cases where reference base is ambiguous
            
            #calls info
            calls=np.fromstring(lineinfo[4], dtype=np.int8) #to ASCII
            
            #qual info
            bq=np.fromstring(lineinfo[5], dtype=np.int8) # base quality, BAQ corrected, ASCII
            mq=np.fromstring(lineinfo[6], dtype=np.int8) # mapping quality, ASCII
            td=np.fromstring(lineinfo[7], dtype=int, sep=',') # distance from tail, comma sep
            
            #find starts of reads ('^' in mpileup)
            startsk=np.where(calls==94)[0]
            for k in startsk:
                calls[k:k+2]=-1
                #remove mapping character, 
                #absolutely required because the next chr could be $
            
            #find ends of reads ('$' in mpileup)
            endsk=np.where(calls==36)[0]
            temp[37]=len(startsk) + len(endsk)
            calls[endsk]=-1
            
            #find indels + calls from reads supporting indels ('+-')
            indelk = np.where((calls==43) | (calls==45))[0]
            for k in indelk:
                forward_looking_index_for_indel_size = 1
                while calls[k + forward_looking_index_for_indel_size] >=48 and calls[k + forward_looking_index_for_indel_size] <58: ## multiple digit indel
                    forward_looking_index_for_indel_size+=1
                if forward_looking_index_for_indel_size > 1:
                    indel_ascii=list(calls[k+1:k+forward_looking_index_for_indel_size])
                    indelsize=int(''.join(map(chr,indel_ascii)))
                    indeld=forward_looking_index_for_indel_size-1
                else:
                    indelsize=int(chr(calls[k+1]))
                    indeld=1
                #record that indel was found in +/- indelregion nearby
                #indexing is slightly different here from matlab version
                if calls[k]==45: #deletion
                    if (position-indelregion-1 >= 0) and (position+indelsize+indelregion-1 < genome_length): # if in middle of contig
                        #must store directly into data as it affects lines earlier and later
                        data[position-indelregion-1:position+indelsize+indelregion-1,39]+=1
                    elif position-indelregion >= 0: # if at end of contig
                        data[position-indelregion-1:,39]+=1
                    else: # if at beginning of contig
                        data[:position+indelsize+indelregion-1,39]+=1
                else: #insertion
                    #insertion isn't indexed on the chromosome, no need for complex stuff
                    if (position-indelregion-1 >= 0) and (position+indelregion-1 < genome_length): # if in middle of contig
                        data[position-indelregion-1:position+indelregion-1,38]+=1
                    elif position-indelregion >= 0: # if at end of contig
                        data[position-indelregion-1:,38]+=1
                    else: # if at beginning of contig
                        data[:position+indelregion-1,38]+=1 # indelsize->indelregion 2022.10.24 Evan and Arolyn

                #remove indel info from counting
                calls[k:(k+1+indeld+indelsize)] = -1 #don't remove base that precedes an indel
            #replace reference matches (.,) with their actual calls
            if ref != None: # when reference allele is not ambiguous
                calls[np.where(calls==46)[0]]=nts_ascii[ref] #'.'
                calls[np.where(calls==44)[0]]=nts_ascii[ref+4] #','
            else: # added 2022.09.23 by Arolyn: for cases where reference allele is ambiguous, confirm there are no ,'s or .'s
                if np.any(calls==46) | np.any(calls==44):
                    print( 'Line from mpileup: ' + line )
                    raise ValueError('Error! Calls at this position allegedly match reference allele even though reference allele was ambiguous.')
            
            #index reads for finding scores
            simplecalls=calls[np.where(calls>0)]
            #simplecalls is a tform of calls where each calls position
            #corresponds to its position in bq, mq, td
            #count how many of each nt and average scores
            for nt in range(8):
                current_nt_indices=np.where(simplecalls == nts_ascii[nt])[0]
                nt_count=len(current_nt_indices)
                if nt_count > 0:
                    temp[nt]=nt_count
                    temp[nt+8]=ghf.round_half_up(np.sum(bq[(current_nt_indices)])/nt_count)-Phred_offset
                    temp[nt+16]=ghf.round_half_up(np.sum(mq[(current_nt_indices)])/nt_count)-33
                    temp[nt+24]=ghf.round_half_up(np.sum(td[(current_nt_indices)])/nt_count)
            
            ## The following section is critical for metagenomic samples 
            # find major and nextmajor allele
            allele_index_summed_calls=[(x,temp[x]+temp[x+4]) for x in range(4)] # add T and t calls, C and c calls, etc. keep index of major allele as first index in tuple
            allele_index_summed_calls.sort(key=lambda x: x[1]) ## sort by number of calls for a given base
            n1 = allele_index_summed_calls[-1][0] # major allele (most calls) (returns actual index of base in nts)
            n2 = allele_index_summed_calls[-2][0] # next most abundant allele (returns actual index of base in nts)
            
            x=np.logical_or(simplecalls==nts_ascii[n1], simplecalls==nts_ascii[n1+4]) # boolean array of reads with major alelle
            y=np.logical_or(simplecalls==nts_ascii[n2], simplecalls==nts_ascii[n2+4]) # boolean array of reads with minor alelle
            
            if sum(temp[0:4]) > min_reads_on_strand and sum(temp[4:8])> min_reads_on_strand and sum(y)*200>sum(x): # only calcualte p values of there are greater than 20 reads on each strand and MAF < .995
                bttests_x=bq[x]-Phred_offset
                bttests_y=bq[y]-Phred_offset # gets base qualities for major/minor alleles
                mttests_x=mq[x]-Phred_offset
                mttests_y=mq[y]-Phred_offset # gets mapping qualities for major/minor alleles
                fttests_x=td[(np.where(simplecalls==nts_ascii[n1])[0])]
                fttests_y=td[(np.where(simplecalls==nts_ascii[n2])[0])] # gets tail distances from fwd reads for major/minor alleles
                rttests_x=td[(np.where(simplecalls==nts_ascii[n1+4])[0])]
                rttests_y=td[(np.where(simplecalls==nts_ascii[n2+4])[0])] # gets tail distances from rev reads for major/minor alleles
                
                ## NOTE: matlab ttest will output result if len(y) == 1, treating it as the hypothesized mean of the normal distribution compared to x 
                bp=ttest_ind(bttests_x,bttests_y)[1]
                mp=ttest_ind(mttests_x,mttests_y)[1]
                fp=ttest_ind(fttests_x,fttests_y)[1]
                rp=ttest_ind(rttests_x,rttests_y)[1]
                # report pval of above t-tests as -log10
                    ## coverting values of 0 and nan for output as -1:
                if bp == 0 or np.isnan(bp):
                    temp[33]=-1
                else: temp[33]=ghf.round_half_up(-log10(bp)) 

                if mp == 0 or np.isnan(mp):
                    temp[34]=-1
                else: temp[34]=ghf.round_half_up(-log10(mp))

                if fp == 0 or np.isnan(fp):
                    temp[35]=-1
                else: temp[35]=ghf.round_half_up(-log10(fp))
                
                if rp == 0 or np.isnan(rp):
                    temp[36]=-1
                else: temp[36]=ghf.round_half_up(-log10(rp))
                
                ## currently not broken but testing against matlab script fails since matlab script *is* broken.
                p=fisher_exact(np.array([[temp[n1], temp[n2]],[temp[n1+4],temp[n2+4]]]), alternative='two-sided')[1] # fisher exact test for strand bias (contingency table = # of reads that are fwd/rev and support major/minor allele)
                if p == 0:
                    temp[32]=-1
                else:
                    temp[32]=ghf.round_half_up(-log10(p))
            
            ## store the data
            #-1 is needed to turn 1-indexed positions to python 0-indexed
            data[position-1,:38]=temp[:38]

    #calc coverage
    coverage=np.sum(data[:,0:8], axis = 1)

    return data, coverage

#%%
if __name__ == "__main__":
    
    SCRIPTS_DIR="scripts"
    sys.path.insert(0, SCRIPTS_DIR)
        
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', dest='input', type=str, help='Path to input pileup',required=True)
    parser.add_argument('-r', dest='ref', type=str, help='Path to reference genome',required=True)
    parser.add_argument('-o', dest='output', type=str, help='Path to output diversity file', required=True)
    parser.add_argument('-c', dest='coverage', type=str, help='Path to coverage file', required=True)
    
    args = parser.parse_args()
    
    diversity_arr, coverage_arr = pileup2diversity(args.input,args.ref)
    
    np.savez_compressed(args.output, data = diversity_arr.astype(int))
    if args.coverage:
        np.savez_compressed(args.coverage, coverage = coverage_arr.astype(int))