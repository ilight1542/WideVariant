#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing helper scripts for GUS
"""
from Bio import SeqIO
from itertools import compress
import os
import numpy as np
import csv
import glob
import subprocess
import gzip
import sys
from math import floor


def get_nt_order():
    order='ATCG'
    return order

def read_samples_CSV(spls,quiet=False):
    """
    Read sample information from a CSV file.

    Parameters:
        spls (str): Path to the CSV file containing sample information.
        quiet (bool, optional): If True, suppresses printing of certain messages. Defaults to False.

    Returns:
        list: A list containing the values from the 'Path', 'Sample', 'FileName', 'Reference', 'Group', and 'Outgroup' columns.

    Raises:
        None

    Description:
        This function reads sample information from a CSV file and returns the values from specific columns.
        The function assumes the CSV file has a specific format and column order.
        The function checks if the first line of the file is a header line and verifies the column order based on it.
        If the header line is missing or the column order is incorrect, warnings or error messages are displayed.
        The function builds lists for each column using the provided header-to-column mapping.
        Finally, the function returns a list containing the values from the specified columns.

    Example:
        sample_info = read_samples_CSV('sample_info.csv', quiet=True)
    """
    smpl_csv_dict = {'Path': [],'Sample': [],'FileName': [],'Reference': [],'Group': [],'Outgroup': []}

    with open(spls, 'r') as fid: 
        for line_id, line in enumerate(fid):
            # check if first line is header. Note: Even when header is not found code continues (w/ warning).
            if line_id == 0:
                if line.startswith('Path,Sample,'):
                    hdr_colnames = line.strip('\n').split(',')
                    hdr_to_col_id = {hdr: id for id, hdr in enumerate(hdr_colnames)}
                    if not quiet:
                        print("Passed CSV header check")
                    continue
                else:
                    Warning("\n\nCSV did NOT pass header check!")
                    if len(line.strip('\n').split(',')) == len(smpl_csv_dict.keys()):
                        Warning('Assumes column order to be: Path,Sample,FileName,Reference,Group,Outgroup\n\n')
                        hdr_to_col_id = {hdr: id for id, hdr in enumerate(smpl_csv_dict.keys())}
                    else:
                        print('\n\nSample csv is not formatted correctly.')
                        print('Snakemake will not start!\n\n')
                        sys.exit(1)
            
            line = line.strip('\n').split(',')
            # build lists
            for key, id in hdr_to_col_id.items():
                smpl_csv_dict[key].append(line[id])
        
    return [smpl_csv_dict[key] for key in smpl_csv_dict.keys()]


def split_samplesCSV(PATH_ls,SAMPLE_ls,FILENAME_ls,REF_Genome_ls,GROUP_ls,OUTGROUP_ls, outdir):
    '''Take info from samples.csv, concat by sample name & save each line as sample_info.csv in data/{sampleID}'''
    
    #Loop through unique samples
    for sample in set(SAMPLE_ls):
        # Concat info for this sample
        sample_info_bool=[s==sample for s in SAMPLE_ls]
        sample_paths=list(compress(PATH_ls,sample_info_bool))
        sample_filenames=list(compress(FILENAME_ls,sample_info_bool))
        sample_references=list(compress(REF_Genome_ls,sample_info_bool))
        sample_groups=list(compress(GROUP_ls,sample_info_bool))
        sample_outgroups=list(compress(OUTGROUP_ls,sample_info_bool))
        
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

def chrpos2index(chrpos,chr_starts):
    '''
    Args:
        chrpos (arr): px2 array of chromsome idx and position on chromosome.
        chr_starts (arr): Vector of chromosome starts (begins at 0).

    Returns:
        p (arr): Vector of position indexes.

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

    return p
    
def p2chrpos(p, ChrStarts):
    '''Convert 1col list of pos to 2col array with chromosome and pos on chromosome

    Args:
        p (array (1 x positions) ): 0-based index of candidate variant positions from Candidate Mutation Table (CMT).
        ChrStarts (array (1 x scaffolds) ): start indices (0-based) of scaffolds in reference genome.

    Returns:
        chrpos (array): DESCRIPTION.

    '''
        
    # get chr and pos-on-chr
    chromo = np.zeros(len(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chromo = chromo + (p >= i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chromo] # [chr-1] -1 due to 0based index
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
    - genome_chr_starts (list): List of starting positions of chromosomes/scaffolds in the genome.
    - scaf_names (list): List of scaffold names in the reference genome.

    Returns:
    - int: Absolute genomic position corresponding to the input chromosome and position.

    Raises:
    - ValueError: If the chromosome_id is not found in the provided scaffold names.
    """
    
    if len(genome_chr_starts) == 1:
        position=int(nt_pos)
    else:
        if chromosome_id not in scaf_names:
            raise ValueError("Scaffold name in pileup file not found in reference")
        position=int(genome_chr_starts[np.where(chromosome_id==scaf_names)]) + int(nt_pos)
    return position

########################################
## Nucleotide str converstion 
########################################

def define_nt_order():
    """
    Define order of nucleotides for converstion of string to integers 
    """
    return 'ATCGatcg'

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
    This function rounds a number to the specified number of decimals using the "round half up" method,
    where .5 values are always rounded up (instead of to even, see 'bankers rounding')

    Parameters:
    - n (float): The number to be rounded.
    - decimals (int, optional): The number of decimals to round to. Default is 0.

    Returns:
    - float: The rounded number.
    """
    multiplier = 10 ** decimals
    return floor(n*multiplier + 0.5) / multiplier  



#TODO: remove if not necessary

# def get_clade_wildcards(cladeID):
#     is_clade = [int(i == cladeID) for i in GROUP_ls]
#     sampleID_clade = list(compress(SAMPLE_ls,is_clade))
#     reference_clade = list(compress(REF_Genome_ls,is_clade))
#     outgroup_clade = list(compress(OUTGROUP_ls,is_clade))
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

