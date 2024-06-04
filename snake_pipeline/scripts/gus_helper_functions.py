#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing helper scripts for GUS
"""
from Bio import SeqIO
from itertools import compress
import os
import numpy as np
import pandas as pd
import csv
import glob
import subprocess
import gzip
import sys
from math import floor


def get_nt_order():
    """
    This function returns a string representing the default order of nucleotides for the entire analysis: 
    Adenine (A), Thymine (T), Cytosine (C), and Guanine (G).

    Returns:
        str: A string 'ATCG' representing the default order of nucleotides.
    """
    order='ATCG'
    return order

def read_samples_CSV(samplescsv_path, expected_columnnames = ['Path','Sample','FileName','Reference','Group', 'Outgroup'],quiet=False):
    """
    Read a CSV file, validate its structure, and return its contents as a dictionary.

    This function performs the following steps:
    1. Checks if the CSV file contains the expected column names (in any order).
    2. Reads the CSV file, with or without headers as necessary.
    3. Validates the CSV file by checking for missing values and duplicated lines.
    4. Removes any rows with missing values or duplicated lines, if found.
    5. Converts the DataFrame to a dictionary and returns it.

    Args:
        samplescsv_path (str): The path to the CSV file to be read.
        expected_columnnames (list of str, optional): The list of expected column names. 
            Defaults to ['Path', 'Sample', 'FileName', 'Reference', 'Group', 'Outgroup'].
        quiet (bool, optional): If True, suppresses print statements for missing values and duplicates. Defaults to False.

    Returns:
        dict: A dictionary representation of the CSV file where keys are column names and values are lists of column data.

    Raises:
        SystemExit: If the CSV file does not contain the expected number of columns.
    """
    # Check if the file contains all expected column names (in any order)
    with open(samplescsv_path, 'r') as file:
        first_line = file.readline().strip().split(',')
    header_present = all([expected_col in first_line for expected_col in expected_columnnames])
    
    if header_present:
        samplescsv = pd.read_csv(samplescsv_path)
    else:
        samplescsv = pd.read_csv(samplescsv_path, header=None)
        if samplescsv.shape[1] != len(expected_columnnames):
            print("samples.csv does not contain the expected number of columns.")
            return sys.exit()
        samplescsv.columns = expected_columnnames
    # Check for missing values & duplicated lines 
    if samplescsv.isna().any().any():
        if not quiet:
            print(f"samples.csv contains missing values. Removing the lines from file:\n{samplescsv[samplescsv.isna().any(axis = 1)].to_string(index=False)}")
        samplescsv = samplescsv[~( samplescsv.isna().any(axis = 1) )]
    if samplescsv.duplicated().any():
        if not quiet:
            print(f"samples.csv contains duplicated lines. Removing one of the duplicated lines from file:\n{samplescsv[samplescsv.duplicated()].to_string(index=False)}")
        samplescsv = samplescsv[~( samplescsv.duplicated() )]
    samplescsv[['Path','Sample','FileName','Reference','Group']] = samplescsv[['Path','Sample','FileName','Reference','Group']].astype(str)
    samplescsv['Outgroup'] = samplescsv['Outgroup'].astype(int)
    if not quiet:
        print(f'{samplescsv.shape[0]} samples are stated in the samples.csv and will be processed.')
    samplescsv_dict = samplescsv.to_dict(orient='list')
    return samplescsv_dict


def split_samplesCSV(samplescsv_dict, outdir):
    """
    This function processes a dictionary containing sample information (read from samples.csv), 
    groups the data by sample name, and saves each group's information into a separate CSV file 
    located in a directory specific to that sample (data/{sampleID}).
    
    Args:
        samplescsv_dict (dict): A dictionary where keys are column names and values are lists of column data.
                                Expected keys are 'Path', 'Sample', 'FileName', 'Reference', 'Group', and 'Outgroup'.
        outdir (str): The directory where the sample-specific subdirectories and CSV files will be saved.
    """
    
    '''Take info from samples.csv, concat by sample name & save each line as sample_info.csv in data/{sampleID}'''
    
    #Loop through unique samples
    for sample in set(samplescsv_dict['Sample']):
        # Concat info for this sample
        sample_info_bool=[s==sample for s in samplescsv_dict['Sample']]
        sample_paths=list(compress(samplescsv_dict['Path'],sample_info_bool))
        sample_filenames=list(compress(samplescsv_dict['FileName'],sample_info_bool))
        sample_references=list(compress(samplescsv_dict['Reference'],sample_info_bool))
        sample_groups=list(compress(samplescsv_dict['Group'],sample_info_bool))
        sample_outgroups=list(compress(samplescsv_dict['Outgroup'],sample_info_bool))
        
        sample_info_csv=list(zip(sample_paths,[sample]*sum(sample_info_bool),sample_filenames,sample_references,sample_groups,sample_outgroups))
        
        path_to_sample_info_csv = f'{outdir}{sample}/sample_info.csv'
        # make data directory for this sample if it doesn't already exist
        if not(os.path.isdir(f'{outdir}{sample}')):
            os.makedirs(f'{outdir}{sample}', exist_ok=True)
        # check to see if this mini csv with sample info already exists
        if os.path.isfile(path_to_sample_info_csv):
            # if so, read file
            old_file = open(path_to_sample_info_csv,'r')
            old_info_read = csv.reader(old_file)
            old_info = list(map(tuple, old_info_read))
            old_file.close()
            
            # check to see if the existing file is consistent with samples.csv
            if not(old_info == sample_info_csv):
                # if not, remove the old file and save sample info in a new file
                # print('Information file must be updated.')
                os.remove(path_to_sample_info_csv)
                with open(path_to_sample_info_csv, "w") as f:
                    writer = csv.writer(f)
                    for row in sample_info_csv:
                        writer.writerow(row)
        else: # if mini csv with sample info does not already exist
            # save sample info in mini csv
            with open(path_to_sample_info_csv, "w") as f:
                writer = csv.writer(f)
                for row in sample_info_csv:
                    writer.writerow(row)


########################################
## Make data links functions
########################################

def makelink(path,sample,filename,output_dir):
    """
    Create symbolic links to FASTQ files for a given sample.

    Parameters:
        path (str): Path to the directory containing the FASTQ files or a file name in batch.
        sample (str): Name of the sample.
        filename (str): Name of the FASTQ file.
        output_dir (str): Output directory where the symbolic links will be created.

    Returns:
        None

    Raises:
        None

    Description:
        This function creates symbolic links to the specified FASTQ files for a given sample in the specified output directory.
        The function first calls the 'findfastqfile' function to locate the FASTQ files based on the provided path, sample, and filename.
        If the FASTQ files are not located using an absolute path or a path relative to the user's home directory, the current working directory is used as the base path.
        The function then prints and executes the commands to create symbolic links using the 'ln -s -T' command.

    Example:
        makelink('/path/to/files', 'sample1', 'file.fastq', '/output/directory')
    """
    [fwd_file, rev_file]=findfastqfile(path,sample, filename)
    if not fwd_file.startswith('/') and not fwd_file.startswith('~/'):
        fwd_file = os.getcwd() + '/' + fwd_file
    if not rev_file.startswith('/') and not rev_file.startswith('~/'):
        rev_file = os.getcwd() + '/' + rev_file
    print(f"ln -s -T {fwd_file} {output_dir}/{sample}/R1.fq.gz")   
    print(f"ln -s -T {rev_file} {output_dir}/{sample}/R2.fq.gz")
    subprocess.run(f"ln -s -T {fwd_file} {output_dir}/{sample}/R1.fq.gz", shell=True)    
    subprocess.run(f"ln -s -T {rev_file} {output_dir}/{sample}/R2.fq.gz", shell=True)  


########################################
## Fastq file functions
########################################

def findfastqfile(dr,ID,filename):
    fwd=[]
    rev=[]
    potentialhits_forward=glob.glob(dr + '/' + filename +'/*1.fastq.gz')
    potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2.fastq.gz')
    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
        fwd=potentialhits_forward[0]
        rev=potentialhits_reverse[0]
    elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0: ## need or statement to further screen path if just one file was found i.e. for *R1_001.fastq.gz *R2_001.fastq.gz
        potentialhits_forward=glob.glob(dr + '/' + filename +'/*1_001.fastq.gz')
        potentialhits_reverse=glob.glob(dr + '/' + filename +'/*2_001.fastq.gz')
        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
            fwd=potentialhits_forward[0]
            rev=potentialhits_reverse[0]
        elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
            potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq.gz')
            potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq.gz')
            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                fwd=potentialhits_forward[0]
                rev=potentialhits_reverse[0]
            elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                potentialhits_forward=glob.glob(dr + '/' + filename +'*1_001.fastq.gz')
                potentialhits_reverse=glob.glob(dr + '/' + filename +'*2_001.fastq.gz')
                if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                    fwd=potentialhits_forward[0]
                    rev=potentialhits_reverse[0]
                elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                    potentialhits_forward=glob.glob(dr + '/' + filename +'*1.fastq')
                    potentialhits_reverse=glob.glob(dr + '/' + filename +'*2.fastq')
                    if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                        fwd=potentialhits_forward[0]
                        rev=potentialhits_reverse[0]
                    elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                        potentialhits_forward=glob.glob(dr + '/' + filename +'*1_001.fastq')
                        potentialhits_reverse=glob.glob(dr + '/' + filename +'*2_001.fastq')
                        if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                            fwd=potentialhits_forward[0]
                            rev=potentialhits_reverse[0]
                    else:
                        foldername=glob.glob(dr + '/' + filename + '*')
                        if foldername and os.path.isdir(foldername[0]):
                            foldername=foldername[0]
                            potentialhits_forward=glob.glob(foldername + '/*' + filename + '*1*.fastq.gz')
                            potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq.gz')
                            if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                                fwd=potentialhits_forward[0]
                                rev=potentialhits_reverse[0]
                            elif len(potentialhits_forward)==0 or len(potentialhits_reverse)==0:
                                print(foldername + '/*' + filename + '*2*.fastq.gz')
                                potentialhits_forward=glob.glob(foldername +  '/*' + filename + '*1*.fastq')
                                potentialhits_reverse=glob.glob(foldername + '/*' + filename + '*2*.fastq')
                                if len(potentialhits_forward)==1 and len(potentialhits_reverse)==1:
                                    fwd=potentialhits_forward[0]
                                    rev=potentialhits_reverse[0]
    if not(fwd) or not(rev):
        raise ValueError('Either no file or more than 1 file found in ' + dr + ' for ' + ID)
    #TODO: add search pattern used to find fastqs for more useful error report
    #TODO: add error differentiation for no file or more than one file  
    
    ##zip fastq files if they aren't already zipped
    if fwd[-3:] != '.gz':
        subprocess.run("gzip " + fwd, shell=True) 
        fwd=fwd+'.gz'
    if rev[-3:] != '.gz':
        subprocess.run("gzip " + rev, shell=True)
        rev=rev+'.gz'
    return [fwd, rev]

  
def cp_append_files(paths,sample,filename,output_dir):
    #When sample is run on multiple lanes with same barcode
    fwd_list=''
    rev_list=''
    for path in paths:
        #Provider name can be either a COMPLETE directory name or a file name in batch(called path in this fx)
        [fwd_file, rev_file]=findfastqfile(path,sample, filename)
        fwd_list=fwd_list+ ' ' + fwd_file
        rev_list=rev_list+ ' ' + rev_file
        print(rev_list)
        print(fwd_list)
    subprocess.run(f"zcat {fwd_list} | gzip > {output_dir}/{sample}/R1.fq.gz", shell=True)
    subprocess.run(f"zcat {rev_list} | gzip > {output_dir}/{sample}/R2.fq.gz", shell=True)


########################################
## Chromosome functions
########################################

def read_fasta(REFGENOME_DIR): 
    '''Reads in fasta file. If directory is given, reads in dir/genome.fasta
    Args:
        REFGENOME_DIR (str): Path to reference genome.

    Returns: SeqIO object for reference genome.
    '''
    fasta_file = glob.glob(REFGENOME_DIR + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOME_DIR + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOME_DIR)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    
    return refgenome

def genomestats(REFGENOMEFOLDER):
    '''Parse genome to extract relevant stats

    Args:
        REFGENOMEFOLDER (str): Directory containing reference genome file.

    Returns:
        ChrStarts (arr): DESCRIPTION.
        Genomelength (arr): DESCRIPTION.
        ScafNames (arr): DESCRIPTION.

    '''

    refgenome = read_fasta(REFGENOMEFOLDER)
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return ChrStarts,Genomelength,ScafNames

def chrpos2index(chrpos,chr_starts,genomelength=np.nan):
    '''
    Args:
        chrpos (arr): 2D array of chromsome ID (1-based) and position on chromosome (0-based).
        chr_starts (arr): Vector of chromosome starts (begins at 0).
        genomelength (float, optional): Total genome length. If specified, checks if any positions exceed this length. Default is np.nan.

    Returns:
        p (arr): Vector of position indexes. (0-based)

    '''
    if (np.size(chrpos,1) != 2) & (np.size(chrpos,0) == 2):
        chrpos=chrpos.T
        print('Reversed orientation of chrpos')
    elif (np.size(chrpos,1) != 2) & (np.size(chrpos,0) != 2):
        print('Wrong input format given. Provide np.array with shape (p, 2)')
        return 
    
    if len(chr_starts) == 1:
        p=chrpos[:,1]
    else:
        p=chr_starts[chrpos[:,0]-1]+chrpos[:,1]
    if not np.isnan(genomelength) and any(p>=genomelength):
        print('Invalid genome positions given with positions larger than genome length')
        return sys.exit()

    return p

def p2chrpos(p, ChrStarts):
    '''Convert 1col list of pos to 2col array with chromosome and pos on chromosome

    Args:
        p (array (1 x positions) ): 0-based index of candidate variant positions from Candidate Mutation Table (CMT).
        ChrStarts (array (1 x scaffolds) ): start indices (0-based) of scaffolds in reference genome.

    Returns:
        chrpos (array): 2D array of chromsome ID (1-based) and position on chromosome (0-based).

    '''
        
    # get chr and pos-on-chr
    chromo = np.ones(np.shape(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chromo = chromo + (p >= i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chromo-1] # [chr-1] -1 due to 0based index
        chrpos = np.column_stack((chromo,positions))
    else:
        chrpos = np.column_stack((chromo,p))
    return chrpos

def convert_chrpos_to_abspos(chromosome_id, nt_pos, genome_chr_starts, scaf_names):  
    """
    convert position on chromosome to absolute position
    Args:
    - chromosome_id (str): Identifier of the chromosome or scaffold of the position which should be converted.
    - nt_pos (int): Position on the chromosome (0-based) which should be converted to absolute position.
    - genome_chr_starts (array): Array of starting positions of chromosomes/scaffolds in the genome.
    - scaf_names (array): Array of scaffold names in the reference genome.

    Returns:
    - int: Absolute genomic position corresponding to the input chromosome and position (0-based).

    Raises:
    - ValueError: If the chromosome_id is not found in the provided scaffold names.
    """
    
    if len(genome_chr_starts) == 1:
        position=int(nt_pos)
    else:
        if chromosome_id not in scaf_names:
            raise ValueError("Scaffold name in pileup file not found in reference")
        index_of_chrom=np.where(scaf_names==chromosome_id)[0]
        position=int(genome_chr_starts[index_of_chrom]) + int(nt_pos)
        if index_of_chrom + 1 < len(scaf_names):
            if position >= genome_chr_starts[index_of_chrom+1]:
                raise ValueError(f"Position {nt_pos} given for chrom {chromosome_id} lies outside the end of the chromosome.")
    return position

########################################
## Nucleotide str converstion 
########################################

def define_nt_order():
    """
    Define order of nucleotides for converstion of string to integers 
    """
    full_nt_order=get_nt_order()+get_nt_order().lower()
    return full_nt_order

def generate_ref_to_int_converstion_dict(nts=define_nt_order()):
    """
    Generate a dictionary for converting nucleotide symbols to integer.

    Parameters:
    - nts (list, optional): List of nucleotide symbols in the desired order. Default is the order defined by `define_nt_order()`.

    Returns:
    - dict: A dictionary where keys are nucleotide symbols, and values are their corresponding integer indices.
    """
    return {nt: i for i, nt in enumerate(nts)}

def convert_ref_to_int(ref_str, nts_dict=generate_ref_to_int_converstion_dict()):
    """
    This function takes a reference allele symbol, e.g. 'A', 'T', 'C', 'G', 'a', 't', 'c', 'g',
    and converts it to its corresponding integer using a provided nucleotide dictionary.

    Parameters:
    - ref_str (str): Reference allele symbol to be converted.
    - nts_dict (dict, optional): Dictionary mapping nucleotide symbols to their integer indices.
                                 Default is the order defined by `generate_ref_to_int_converstion_dict()`.

    Returns:
    - int or None: The integer index corresponding to the reference allele symbol.
                  Returns None for cases where the reference base is ambiguous.
    """
    if ref_str in nts_dict.keys():
        ref=nts_dict[ref_str] # convert to 0123
        if ref >= 4:
            ref = ref - 4
    else:
        ref=None # for cases where reference base is ambiguous
    return ref


########################################
## Miscellaneous functions
########################################

def round_half_up(n, decimals=0):
    """
    Round a number to the nearest integer using half-up rounding.

    This function rounds the input number `n` to the nearest integer using half-up rounding. 
    Half-up rounding means that if the fraction part is exactly 0.5, it rounds up to the next integer.

    Args:
        n (float): The number to be rounded.
        decimals (int, optional): The number of decimal places to round to (default is 0).

    Returns:
        float: The rounded number.

    """
    multiplier = 10 ** decimals
    return floor(n*multiplier + 0.5) / multiplier   



#TODO: remove if not necessary

# def get_clade_wildcards(cladeID):
#     is_clade = [int(i == cladeID) for i in samplescsv_dict['Group']]
#     sampleID_clade = list(compress(samplescsv_dict['Sample'],is_clade))
#     reference_clade = list(compress(samplescsv_dict['Reference'],is_clade))
#     outgroup_clade = list(compress(samplescsv_dict['Outgroup'],is_clade))
#     return sampleID_clade,reference_clade,outgroup_clade
    
# def get_sampleID_names(wildcards):  
#     sampleID_clade,_,_ = get_clade_wildcards(wildcards.cladeID)
#     return sampleID_clade

# def get_outgroup_bool(wildcards):  
#     _,_,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     return outgroup_clade

# def get_positions_prep(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     mat_positions_prep=expand("Case/temp/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.npz",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return mat_positions_prep

# def get_diversity(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     diversity_mat = expand("Mapping/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.npz",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return diversity_mat   

# def get_quals(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     quals_mat = expand("Mapping/quals/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.npz",zip,sampleID=sampleID_clade, reference=reference_clade, outgroup=outgroup_clade)
#     return quals_mat 

# def get_ref_genome(wildcards):
#     sampleID_clade,reference_clade,outgroup_clade = get_clade_wildcards(wildcards.cladeID)
#     return reference_clade

