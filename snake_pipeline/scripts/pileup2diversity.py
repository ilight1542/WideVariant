#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 22:14:06 2022

@author: evanqu
"""

import numpy as np
import sys
import argparse
import gus_helper_functions as ghf
from scipy.stats import ttest_ind, fisher_exact
from math import log10

#%% Version history
#2022.02.08: Evan: Direct translation from pileup_to_diversity_matrix_snakemake.m
#2022.10.18, Arolyn: Now works when reference genome has lowercase letters or ambiguous letters
#2022.10.23, Arolyn: Updated comments on 40 statistics to have python indexing (0-39) as opposed to matlab indexing (1-40)

#%%Some notes

# This function saves the following 40 statistics for each position on the
# reference genome for the sample being analyzed:
# [A T C G a t c g bqA bqT ... bqg mqA .... mqg  tdA .... tdg Ps Pb Pm Pftd Prtd E I D]
# List of statistics by index:
# [0-3] A is the number of forward reads supporting A
# [4-7] a is the number of reverse reads supporting A
# [8-15] bq is the average phred qualities of all A's
# [16-23] bm is the average mapping qualities of all A's
# [24-31] td is the average tail distance of all A's
# [32] Ps is the p value for strand bias (fishers test)
# [33] Pb is the p value for the base qualities being the same for the two
# different types of calls (1st major, 2nd major nt, either strand) (ttest)
# [34] Pm is the p value for the mapping qualities being the same for the two
# different types of calls (ttest)
# [35] Pftd is the p value for the tail distantces on the forward strand
# being the same for the two different types of calls (ttest)
# [36] Pftd is the p value for the tail distantces on the reverse strand
# being the same for the two  different types of calls (ttest)
# [37] E is number of calls at start or ends of a read
# [38] I is number of reads supporting insertions in the +/- (indel_region) bp region
# [39] D is number of reads supporting deletions in the +/- (indel_region) bp region

# ChrStarts: is an array holding the indices in the position dimension
# corresponding to the start of a new chromsome.

#%%

def generate_nts_ascii():
    """
    Generate ASCII representations of nucleotides.

    This function generates ASCII representations of nucleotides by combining uppercase and lowercase nucleotides from the get_nt_order function in the ghf module.

    Returns:
        list: A list of ASCII representations of nucleotides.

    """
    nts=ghf.get_nt_order()+ghf.get_nt_order().lower()
    nts_ascii=[ord(nt) for nt in nts]
    return nts_ascii

def clean_calls_from_start_and_end(calls):
    """
    Clean the calls array by modifying starts and ends of read characters to -1.

    Args:
        calls (ndarray): Array of called bases.

    Returns:
        ndarray: Cleaned array of called bases.

    """
    #find starts of reads ('^' in mpileup)
    startsk=np.where(calls==94)[0]
    for k in startsk:
        calls[k:k+2]=-1
        #remove mapping character, 
        #absolutely required because the next chr could be $
    
    #find ends of reads ('$' in mpileup)
    endsk=np.where(calls==36)[0]
    num_reads_start_end=len(startsk) + len(endsk)
    calls[endsk]=-1
    return calls,num_reads_start_end

def get_indel_size(calls,k):
    """
    Find the total size of a call for a given indel, and the width of the declaration of the indel size from calls.

    Args:
        calls (ndarray): Array of called bases for a given position and sample
        k: index of calls which supports an indel 
            info on indel call follow directly after in format: +11ATTCATTCATT 
                + = insertion beginning directly after this position
                11 = of 11bp
                ATTCATTCATT = this read has an insertion of ATTCATTCATT 
        
    Returns:
        Tuple[int, int]: Tuple containing the size of the indel (defined as 11 in the example above),
            ,and the number of chars in calls used to define the size (which would be 2 since 11 takes up 2 string positions in the pileup calls entry)
    """
    forward_looking_index_for_indel_size = 1
    while calls[k + forward_looking_index_for_indel_size] >=48 and calls[k + forward_looking_index_for_indel_size] <58: ## multiple digit indel
        forward_looking_index_for_indel_size+=1
    if forward_looking_index_for_indel_size > 1:
        indel_ascii=list(calls[k+1:k+forward_looking_index_for_indel_size])
        indelsize=int(''.join(map(chr,indel_ascii)))
        indeld=forward_looking_index_for_indel_size-1 # distance ahead of current call index where indel base call info is stored
    else:
        indelsize=int(chr(calls[k+1]))
        indeld=1
    return indelsize,indeld

def parse_indels_into_data(calls,position,data,indel_region,genome_length):
    """
    Find indels and calls from reads supporting indels.

    Args:
        calls (ndarray): Array of called bases for a given position and sample.
        data (ndarray): Data structure for recording indel information.
        position (int): Current position in the genome.
        indel_region (int): Size of the indel region.
        genome_length (int): Length of the genome.

    Returns:
        Tuple[ndarray, ndarray]: Tuple containing the updated calls array and data structure.

    """
    #find reads supporting indels in calls ('+-')
    indelk = np.where((calls==43) | (calls==45))[0]
    for k in indelk:
        # get number of bases and indexing info for calls of indel
        indelsize,indeld=get_indel_size(calls,k)
        # Mark indel region bases before and after indel
        # record that indel support was found in +/- indel_region nearby
        indel_region_start=position-indel_region+1 # indel starts after current position, so pos 100, with indel support between 100 and 101, gets 98,99 and 100 updated
        del_region_end=position+indelsize+indel_region+1 # indel starts after current position, so pos 100, with indel support between 100 and 101, 101,102,103 updated (plus size of del)
        ins_region_end=position+indel_region+1 # indel starts after current position, so pos 100, with indel support between 100 and 101, 101,102,103 updated (plus size of del)
        if calls[k]==45: #deletion
            if (indel_region_start >= 0) and (del_region_end < genome_length): # if in middle of contig
                #must store directly into data as it affects lines earlier and later
                data[indel_region_start:del_region_end,39]+=1
            elif position-indel_region >= 0: # if at end of contig
                data[indel_region_start:,39]+=1
            else: # if at beginning of contig
                data[:del_region_end,39]+=1
        else: #insertion
            #insertion isn't indexed on the chromosome, no need for complex stuff
            if (indel_region_start >= 0) and (ins_region_end < genome_length): # if in middle of contig
                data[indel_region_start:ins_region_end,38]+=1
            elif position-indel_region >= 0: # if at end of contig
                data[indel_region_start:,38]+=1
            else: # if at beginning of contig
                data[:ins_region_end,38]+=1 # indelsize->indel_region 2022.10.24 Evan and Arolyn

        #mask indel info from calls
        calls[k:(k+1+indeld+indelsize)] = -1 #don't remove base that precedes an indel
    return calls,data

def parse_calls_into_simplecalls(calls,ref_idx,nts_ascii,line,line_id):
    """

    Parse the calls into simple calls and replace reference matches.

    Args:
        calls (ndarray): Array of called bases.
        ref_idx (int or None): Reference allele index (in nts order from ghf.get_nt_order() ). None if the reference allele is ambiguous.
        nts_ascii (list): List of ASCII representations of nucleotides.
        line (str): line from mpileup file
        line_id (int): 0-based line number from mpileup file

    Returns:
        ndarray: Array of simple calls with reference matches replaced and non-call info removed (eg indel or start/end of read marker)

    Raises:
        ValueError: If reference allele is ambiguous and there are reference matches in the calls.
    
    """
    #simplecalls is a tform of calls where each calls position
    #corresponds to its position in bq, mq, td
    #count how many of each nt and average scores
    #replace reference matches (.,) with their actual calls

    if ref_idx != None: # when reference allele is not ambiguous
        calls[np.where(calls==46)[0]]=nts_ascii[ref_idx] #'.' // fwd read 
        calls[np.where(calls==44)[0]]=nts_ascii[ref_idx+4] #',' // rev read
    else: # added 2022.09.23 by Arolyn: for cases where reference allele is ambiguous, confirm there are no ,'s or .'s
        if np.any(calls==46) | np.any(calls==44):
            raise ValueError(f'Error! Calls at this position allegedly match reference allele even though reference allele was ambiguous. Line {str(line_id+1)} from mpileup: {line}')
    
    #index reads for finding scores
    simplecalls=calls[np.where(calls>0)]
    return simplecalls

def run_statistical_tests(simplecalls,temp,bq,mq,td,nts_ascii,major_nt_idx,read_w_major_nt_bool,minor_nt_idx,read_w_minor_nt_bool,Phred_offset=33):
    """
    Run statistical tests on sequencing data for metagenomic samples. (only if MAF <0.995 and each strand has > min_reads_on_strand)

    Args:
        simplecalls (ndarray): Array of called bases for each read.
        temp (ndarray): Array to store the results of the statistical tests.
        bq (ndarray): Base qualities for each read.
        mq (ndarray): Mapping qualities for each read.
        td (ndarray): Tail distances from forward and reverse reads.
        nts_ascii (list): List of ASCII representations of nucleotides.
        Phred_offset (int): Phred quality score offset.

    Returns:
        ndarray: Updated temp array with statistical test results on indicies 33-36 (inclusive) added (-1 means that the data gave a p-val of 0).

    """
    ## The following section is critical for metagenomic samples 
    # find major and nextmajor allele

    bttests_major_nt,bttests_minor_nt = bq[read_w_major_nt_bool]-Phred_offset, bq[read_w_minor_nt_bool]-Phred_offset # gets base qualities for major/minor alleles
    mttests_major_nt,mttests_minor_nt = mq[read_w_major_nt_bool]-Phred_offset, mq[read_w_minor_nt_bool]-Phred_offset # gets mapping qualities for major/minor alleles
    fttests_major_nt,fttests_minor_nt = td[(np.where(simplecalls==nts_ascii[major_nt_idx])[0])], td[(np.where(simplecalls==nts_ascii[minor_nt_idx])[0])] # gets tail distances from fwd reads for major/minor alleles
    rttests_major_nt,rttests_minor_nt = td[(np.where(simplecalls==nts_ascii[major_nt_idx+4])[0])], td[(np.where(simplecalls==nts_ascii[minor_nt_idx+4])[0])] # gets tail distances from rev reads for major/minor alleles
    
    ## NOTE: matlab ttest will output result if len(y) == 1, treating it as the hypothesized mean of the normal distribution compared to x 
    bq_pval=ttest_ind(bttests_major_nt,bttests_minor_nt)[1]
    mq_pval=ttest_ind(mttests_major_nt,mttests_minor_nt)[1]
    forward_pval=ttest_ind(fttests_major_nt,fttests_minor_nt)[1]
    reverse_pval=ttest_ind(rttests_major_nt,rttests_minor_nt)[1]
    # record pval of above t-tests as -log10
    # coverting values of 0 and nan for output as -1 (test failed):
    if bq_pval == 0 or np.isnan(bq_pval):
        temp[33]=-1
    else: temp[33]=ghf.round_half_up(-log10(bq_pval)) 

    if mq_pval == 0 or np.isnan(mq_pval):
        temp[34]=-1
    else: temp[34]=ghf.round_half_up(-log10(mq_pval))

    if forward_pval == 0 or np.isnan(forward_pval):
        temp[35]=-1
    else: temp[35]=ghf.round_half_up(-log10(forward_pval))
    
    if reverse_pval == 0 or np.isnan(reverse_pval):
        temp[36]=-1
    else: temp[36]=ghf.round_half_up(-log10(reverse_pval))
    
    ## currently not broken but testing against matlab script fails since matlab script *is* broken.
    p=fisher_exact(np.array([[temp[major_nt_idx], temp[minor_nt_idx]],[temp[major_nt_idx+4],temp[minor_nt_idx+4]]]), alternative='two-sided')[1] # fisher exact test for strand bias (contingency table = # of reads that are fwd/rev and support major/minor allele)
    if p == 0:
        temp[32]=-1
    else: temp[32]=ghf.round_half_up(-log10(p))
    return temp

def parse_simplecalls_into_temp_data(simplecalls,temp,bq,mq,td,nts_ascii,min_reads_on_strand,max_MAF,Phred_offset):
    """
    Parse simplecalls into temp data structure for recording read support, basequality values, mapping quality values, and tail distance values for this position for each read possible (ATCG and atcg)

    Args:
        simplecalls (ndarray): Array of called bases for each read.
        temp (ndarray): Array to store the results of the statistical tests.
        bq (ndarray): Base qualities for each read.
        mq (ndarray): Mapping qualities for each read.
        td (ndarray): Tail distances from forward and reverse reads.
        nts_ascii (list): List of ASCII representations of nucleotides.
        min_reads_on_strand (int): Minimum number of reads required on a strand.
        Phred_offset (int): Phred quality score offset.

    Returns:
        ndarray: Updated temp array with values recorded on indicies 0-31 (inclusive) added.

    """

    for nt in range(8):
        current_nt_indices=np.where(simplecalls == nts_ascii[nt])[0]
        nt_count=len(current_nt_indices)
        if nt_count > 0:
            temp[nt]=nt_count
            temp[nt+8]=ghf.round_half_up(np.sum(bq[(current_nt_indices)])/nt_count)-Phred_offset
            temp[nt+16]=ghf.round_half_up(np.sum(mq[(current_nt_indices)])/nt_count)-Phred_offset
            temp[nt+24]=ghf.round_half_up(np.sum(td[(current_nt_indices)])/nt_count)
    
    allele_index_summed_calls=[(x,temp[x]+temp[x+4]) for x in range(4)] # add T and t calls, C and c calls, etc. keep index of major allele as first index in tuple
    allele_index_summed_calls.sort(key=lambda x: x[1]) ## sort by number of calls for a given base
    major_nt_idx = allele_index_summed_calls[-1][0] # major allele (most calls) (returns actual index of base in nts)
    minor_nt_idx = allele_index_summed_calls[-2][0] # next most abundant allele (returns actual index of base in nts)
    
    read_w_major_nt_bool=np.logical_or(simplecalls==nts_ascii[major_nt_idx], simplecalls==nts_ascii[major_nt_idx+4]) # boolean array of reads with major allele
    read_w_minor_nt_bool=np.logical_or(simplecalls==nts_ascii[minor_nt_idx], simplecalls==nts_ascii[minor_nt_idx+4]) # boolean array of reads with minor allele
    
    num_reads_major_nt = sum(read_w_major_nt_bool)
    num_reads_minor_nt = sum(read_w_minor_nt_bool)

    if (sum(temp[0:4]) > min_reads_on_strand) and ( sum(temp[4:8]) > min_reads_on_strand ) and (num_reads_minor_nt > (num_reads_major_nt*(1-max_MAF))): # only calcualte p values of there are greater than min reads on each strand and MAF < max_MAF
        temp=run_statistical_tests(simplecalls,temp,bq,mq,td,min_reads_on_strand,major_nt_idx,read_w_major_nt_bool,minor_nt_idx,read_w_minor_nt_bool,Phred_offset)
    return temp

def parse_entry_in_mpileupup(line,line_id,data,chr_starts,genome_length,scaf_names,nts_ascii,indel_region,min_reads_on_strand,maxMAF,Phred_offset,num_fields_diversity_arr):
    lineinfo = line.strip().split('\t')
    
    #holds info for each position before storing in data
    temp = np.zeros((num_fields_diversity_arr),dtype=int)
    
    chromo = lineinfo[0]
    chr_idx = np.where(chromo == scaf_names)[0][0] +1 ## +1 as chrpos2index uses 1-indexed chromosomes
    position_on_chr = int(lineinfo[1]) #1-indexed
    chr_pos_array = np.array( ([chr_idx, position_on_chr], ) ) ## 2D-np array required!
    position = ghf.chrpos2index(chr_pos_array, chr_starts)

    #ref allele
    ref_str = lineinfo[2]; # reference allele from pileup (usually A T C or G but sometimes a different symbol if nucleotide is ambiguous)
    nts=ghf.get_nt_order()
    nts_dict={nt: i for i, nt in enumerate(nts)}
    if ref_str in nts_dict.keys():
        ref_idx=nts_dict[ref_str] # convert to 0123
        if ref_idx >= 4:
            ref_idx = ref_idx - 4
    else:
        ref_idx=None # for cases where reference base is ambiguous
    
    #calls info
    calls=np.frombuffer(lineinfo[4].encode('utf-8'), dtype=np.int8).copy() #to ASCII // copy to make line permutable again 

    #qual info
    bq=np.frombuffer(lineinfo[5].encode('utf-8'), dtype=np.int8).copy() # base quality, BAQ corrected, ASCII // copy to make line permutable again 
    mq=np.frombuffer(lineinfo[6].encode('utf-8'), dtype=np.int8).copy() # mapping quality, ASCII // copy to make line permutable again 
    td=np.array([num for num in lineinfo[7].split(',') if num.isdigit()], dtype = np.int8)
    
    # find start and end reads supporting call, mask for downstream processing
    calls,temp[37] = clean_calls_from_start_and_end(calls)
    
    # parse indels
    calls,data = parse_indels_into_data(calls,position,data,indel_region,genome_length)

    simplecalls = parse_calls_into_simplecalls(calls,ref_idx,nts_ascii,line,line_id)
    temp = parse_simplecalls_into_temp_data(simplecalls,temp,bq,mq,td,nts_ascii,min_reads_on_strand,maxMAF,Phred_offset)
    
    ## store the data
    #-1 is needed to turn 1-indexed positions to python 0-indexed
    data[position-1,:38]=temp[:38]
    return data

def mpileup_to_diversity(input_pileup,path_to_ref,min_reads_on_strand,indel_region,maxMAF,Phred_offset,num_fields_diversity_arr):
    """
    Convert an mpileup file to a diversity data structure.

    Args:
        input_pileup (str): Path to the input mpileup file.
        path_to_ref (str): Path to the reference genome file.
        min_reads_on_strand (int): Minimum number of reads on a strand in order to run statistical tests for detecting bias in heterozygous called sites
        indel_region (int): Number of bases upstream/downstream of an indel to mark as having an indel-supporting read nearby
        Phred_offset (int): Phred quality score offset.
        num_fields_diversity_arr (int): Number of fields in the diversity data structure.

    Returns:
        ndarray: A data structure containing diversity information.

    """
    #####
    [chr_starts,genome_length,scaf_names] = ghf.genomestats(path_to_ref)
    ## generate list of nts as ascii_idx
    nts_ascii = generate_nts_ascii()
    
    #init
    data = np.zeros((genome_length,num_fields_diversity_arr),dtype=int) #format for each base in genome: [ #A, #T, #C, #G, #a, #t, #c, #g,...other statistics...]

    #Read in mpileup file
    print(f"Reading input file: {input_pileup}")
    with open(input_pileup) as mpileup:
        for line_id, line in enumerate(mpileup):
            data = parse_entry_in_mpileupup(line,line_id,data,chr_starts,genome_length,scaf_names,nts_ascii,indel_region,min_reads_on_strand,maxMAF,Phred_offset,num_fields_diversity_arr)
    return data

def main(input_pileup, path_to_ref,phred_offset=33,num_fields_diversity_arr=40,min_reads_on_strand=20,maxMAF=0.995,indel_region=3):
    """Grabs relevant allele info from mpileupfile and stores as an array 

    Args:
        input_pileup (str): Path to input pileup file.
        path_to_ref (str): Path to reference genome file
        phred_offset (int): Phred quality score offset.
        min_reads_on_strand: default=20: Minimum number of reads on a strand in order to run statistical tests for detecting bias in heterozygous called sites
        indel_region: default=3: Number of bases upstream/downstream of an indel to mark as having an indel-supporting read nearby
    """
    
    data = mpileup_to_diversity(input_pileup,path_to_ref,min_reads_on_strand,indel_region,maxMAF,phred_offset,num_fields_diversity_arr)

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
    parser.add_argument('-p', dest='phred_offset', type=int, help='Phred offset of sequencing data', required=False, default=33)
    
    args = parser.parse_args()
    
    diversity_arr, coverage_arr = main(args.input,args.ref,args.phred_offset)
    
    np.savez_compressed(args.output, data = diversity_arr.astype(int))
    if args.coverage:
        np.savez_compressed(args.coverage, coverage = coverage_arr.astype(int))
