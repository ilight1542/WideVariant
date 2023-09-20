# -*- coding: utf-8 -*-
"""
---Gathers everything together for candidate_mutation_table---
NOTE: Still reads in many *.mat files etc. Further purging of matlab necessary!
Output:
# path_candidate_mutation_table: where to write
# candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat
---
# Inputs (changed to argparse usage):
     path_to_p_file: where to find all_positions.mat
     path_to_sample_names_file: where to find text file with sample names
         (space delimited)
     path_to_outgroup_boolean_file: where to find text file with outgroup
         booleans (space delimited, 1=outgroup, 0=not)
    path_to_list_of_quals_files: where to find text file with list of
       quals.mat files for each sample (space delimited)
     path_to_list_of_diversity_files: where to find text file with list of
       diversity.mat files for each sample (space delimited)
# Output:
     path_candidate_mutation_table: where to write
     candidate_mutation_table.mat, ex. results/candidate_mutation_table.mat
# Note: All paths should be relative to pwd!
## Version history
     This is adapted from TDL's build_mutation_table_master_smaller_file_size_backup.m
  #   Arolyn, 2018.12.19: This script was written as part of the transition to snakemake.
          It performs the part of the case step that gathers data for
         Quals and counts and saves candidate_mutation_table.mat
  #   Arolyn, 2019.02.12: Added another matlab variable that stores indel
          statistics called 'indel_counter'.
  #   Tami, 2019.12.12: Converted into python and also added ability save coverage data
  #   Felix: 2020.01-04: Continous Debugged and adapted script for streamlined Snakemake implementation.
  #                      Added argparse for proper argument parsing and optional coverage matrix build.
  #   Evan: 2022.02.07: Changed to work with fully python version of lab pipeline
  #   Arolyn, 2022.10.17:
        Fixed coverage matrix generation
        Changed inputs for coverage matrices from booleans to output filenames
        Changed coverage matrices to numpy arrays instead of scipy sparse matrices (because the matrices shouldn't be sparse)
  #   Arolyn, 2022.10.23:
        Fixed bug where indices of indel statistics were not correct (line 173): 8:10 -> 38:40
"""

# %%Testing
# os.chdir("/Users/evanqu/Dropbox (MIT)/Lieberman Lab/Personal lab notebooks/Evan/1-Projects/strainslicer/dev/build_cmt_py")
# path_to_p_file="6-case-temp/allpositions.pickle"
# path_to_sample_names_file="6-case-temp/string_sampleID_names.txt"
# path_to_outgroup_boolean_file="6-case-temp/positions/string_outgroup_bool.txt"
# path_to_list_of_quals_files="6-case-temp/string_qual_mat.txt"
# path_to_list_of_diversity_files="6-case-temp/string_diversity_mat.txt"
# path_to_candidate_mutation_table="test_candidate_mutation_table.pickle.gz"


# main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files, path_to_list_of_diversity_files, path_to_candidate_mutation_table, cov_mat_raw, cov_mat_norm ):

# %%
''' load libraries '''
import numpy as np
import pickle
import scipy.io as sio
import os
import sys, argparse
import gzip
from scipy import sparse
from pathlib import Path


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
                            Gathers everything together for candidate_mutation_table.
                            Optional: Builds coverage matrix (optional w/ double-standardized matrix)
                               ''',
                                 epilog="Questions or comments? --> fkey@mit.edu")
parser.add_argument("-p", dest="allpositions", help="All positions p file (.npz)", required=True, action='store')
parser.add_argument("-s", dest="sampleNames", help="File with sample names", required=True, action='store')
parser.add_argument("-g", dest="outgroupBool", help="String outgroup bool", required=True, action='store')
parser.add_argument("-q", dest="qualfiles", help="String qual matrix paths", required=True, action='store')
parser.add_argument("-d", dest="divfiles", help="String diversity paths", required=True, action='store')
parser.add_argument("-o", dest="candidate_mutation_table",
                    help="Output candidate mutation table. Py pickle structure (*.pickle.gz)", required=True,
                    action='store')
# parser.add_argument("-c", dest="get_cov", help="Set flag to build raw coverage matrix as sparse csr gzip numpy object (dirname+cov_raw_sparsecsr_mat.npz)",action="store_true", default=False)
# parser.add_argument("-n", dest="get_dbl_norm_cov", help="Set flag to build double normalized coverage matrix as sparse csr gzip numpy object (dirname+cov_norm_sparsecsr_mat.npz)",action="store_true", default=False)
parser.add_argument("-c", dest="cov_mat_raw", help="Output raw coverage matrix as sparse csr gzip numpy object (*.npz)",
                    action='store', default='')
parser.add_argument("-n", dest="cov_mat_norm",
                    help="Output double normalized coverage matrix as sparse csr gzip numpy object (*.npz)",
                    action='store', default='')
parser.add_argument("-t", dest="dim", help="Specify the number of statistics (default 8)", type=int, default=8)
args = parser.parse_args()

# %%
'''Functions'''


def main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files,
         path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_to_cov_mat_raw, path_to_cov_mat_norm,
         flag_cov_raw, flag_cov_norm, dim):
    pwd = os.getcwd()

    # p: positions on genome that are candidate SNPs
    print('Processing candidate SNP positions...')
    infile=np.load(path_to_p_file) # from previous step, should include variable called p
    p=infile['p'].flatten()
    print('Total number of positions: ' + str(len(p)))

    # SampleNames: list of names of all samples
    print('Processing sample names...')
    fname = pwd + '/' + path_to_sample_names_file
    with open(fname, 'r') as f:
        SampleNames = f.read().splitlines()
    numSamples = len(SampleNames)  # save number of samples
    print('Total number of samples: ' + str(numSamples))

    ## in_outgroup: booleans for whether or not each sample is in the outgroup
    print('Processing outgroup booleans...')

    fname = pwd + '/' + path_to_outgroup_boolean_file
    with open(fname, 'r') as f:
        in_outgroup_str = f.read().splitlines()
    
    in_outgroup = np.asarray([s == '1' for s in in_outgroup_str], dtype=bool).reshape(1, len(in_outgroup_str))

    ## Quals: quality score (relating to sample purity) at each position for all samples
    print('Gathering quality scores at each candidate position...')
    # Import list of directories for where to quals for each sample
    fname = pwd + '/' + path_to_list_of_quals_files
    with open(fname, 'r') as f:
        paths_to_quals_files = f.read().splitlines()
    # Make Quals
    Quals = np.zeros((len(p), numSamples), dtype='int')  # initialize
    for i in range(numSamples):
        print('Loading quals matrix for sample: ' + str(i))
        print('Filename: ' + paths_to_quals_files[i])
        infile=np.load(paths_to_quals_files[i]) # from previous step, should include variable called p
        quals=infile['quals'].flatten()
        Quals[:, i] = quals[p - 1]  # -1 convert position to index

    ## counts: counts for each base from forward and reverse reads at each candidate position for all samples
    print('Gathering counts data at each candidate position...\n')

    # Import list of directories for where to diversity file for each sample
    fname = pwd + '/' + path_to_list_of_diversity_files
    with open(fname, 'r') as f:
        paths_to_diversity_files = f.read().splitlines()
    # Load in first diversity to get some stats
    tempfile=np.load(paths_to_diversity_files[1]) 
    data=tempfile['data']
    size = np.shape(data)
    GenomeLength = size[0]

    # Make counts and coverage at the same time
    counts = np.zeros((dim, len(p), numSamples), dtype='uint')  # initialize
    all_coverage_per_bp = np.zeros((GenomeLength, numSamples), dtype='uint')
    indel_counter = np.zeros((2, len(p), numSamples), dtype='uint')


    for i in range(numSamples):
        print('Loading counts matrix for sample: ' + str(i))
        print('Filename: ' + paths_to_diversity_files[i])
        infile=np.load(paths_to_diversity_files[i]) 
        data=infile['data']
        counts[:, :, i] = data[p - 1, 0:dim].T  # -1 convert position to index

        if flag_cov_raw:
            np.sum(data[:, 0:dim], axis=1, out=all_coverage_per_bp[:, i])

        indel_counter[:, :, i] = data[p - 1, 38:40].T  # Num reads supporting indels and reads supporting deletions
        # -1 convert position to index


    # Normalize coverage by sample and then position; ignore /0 ; turn resulting inf to 0

    if flag_cov_norm:
        with np.errstate(divide='ignore', invalid='ignore'):
            # 1st normalization
            array_cov_norm = (all_coverage_per_bp - np.mean(all_coverage_per_bp, axis=1, keepdims=True)) / np.std(
                all_coverage_per_bp, axis=1,
                keepdims=True)  # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
            array_cov_norm[~np.isfinite(array_cov_norm)] = 0

            # 2nd normalization
            array_cov_norm = (array_cov_norm - np.mean(array_cov_norm, axis=0, keepdims=True)) / np.std(array_cov_norm,
                                                                                                        axis=0,
                                                                                                        keepdims=True)  # ,keepdims=True maintains 2D array (second dim == 1), necessary for braodcasting
            array_cov_norm[~np.isfinite(array_cov_norm)] = 0

            # Scale and convert to int to save space
            array_cov_norm_scaled = (np.round(array_cov_norm, 3) * 1000).astype('int64')

    ## get coverage stats
    ## NOTE: Indexes 0-10: sites with coverage 1, 2, 3,... 
    ## NOTE: Index 11: median coverage across all sites
    ## NOTE: Index 12: mean coverage across all sites
    ## NOTE: Index 13: standard deviation of coverage across all sites
    coverage_stats=np.zeros((numSamples,14),dtype='uint') 
    for index in range(numSamples):
        unique,counts_for_dict=np.unique(all_coverage_per_bp[:, index], return_counts=True)
        counts_dict=dict(zip(unique,counts_for_dict))
        ## find mean and median, create bins 0-10 four outputting
        counts_for_output=[0]*11
        total_covg=0
        total_pos=0
        for val in counts_dict:
            if val <= 10:
                counts_for_output[val]=counts_dict[val]
            else:
                counts_for_output[10]+=counts_dict[val]
            total_pos+=counts_dict[val]
            total_covg+=val*counts_dict[val]
        mean_covg=total_covg/total_pos
        median_val=total_pos/2
        count_to_median=0
        index_for_median=0
        while count_to_median < median_val:
            count_to_median+=counts_for_dict[index_for_median]
            index_for_median+=1
        median=index_for_median-1
        ## save into coverage_stats
        coverage_stats[index,0:11]=counts_for_output
        coverage_stats[index,11]=median
        coverage_stats[index,12]=mean_covg
        coverage_stats[index,13]=np.std(all_coverage_per_bp[index])

    # Reshape & save matrices

    Quals = Quals.transpose() #quals: num_samples x num_pos
    p = p.transpose() #p: num_pos
    counts = counts.swapaxes(0,2) #counts: num_samples x num_pos x 8
    in_outgroup = in_outgroup.flatten() #in_outgroup: num_samples 
    SampleNames = SampleNames #sampleNames: num_samples
    indel_counter = indel_counter.swapaxes(0,2) #indel_counter: num_samples x num_pos x 2

    outdir = os.path.dirname(path_to_candidate_mutation_table)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if flag_cov_raw:
        print("Saving " + path_to_cov_mat_raw)
        all_coverage_per_bp = all_coverage_per_bp.transpose() #all_coverage_per_bp: num_samples x num_pos
        with open(path_to_cov_mat_raw, 'wb') as f:
            np.savez_compressed(path_to_cov_mat_raw, all_coverage_per_bp=all_coverage_per_bp)

    if flag_cov_norm:
        print("Saving " + path_to_cov_mat_norm)
        array_cov_norm_scaled = array_cov_norm_scaled.transpose() #array_cov_norm_scaled: num_samples x num_pos
        with open(path_to_cov_mat_norm, 'wb') as f:
            np.savez_compressed(path_to_cov_mat_norm, array_cov_norm_scaled=array_cov_norm_scaled)

    CMT = {'sample_names': SampleNames,
           'p': p,
           'counts': counts,
           'quals': Quals,
           'in_outgroup': in_outgroup,
           'indel_counter': indel_counter,
           'coverage_stats': coverage_stats}

    file_path = path_to_candidate_mutation_table
    print("Saving " + path_to_candidate_mutation_table)
    np.savez_compressed(file_path, **CMT)

    print('DONE')


# %%
if __name__ == "__main__":
    path_to_p_file = args.allpositions
    path_to_sample_names_file = args.sampleNames
    path_to_outgroup_boolean_file = args.outgroupBool
    path_to_list_of_quals_files = args.qualfiles
    path_to_list_of_diversity_files = args.divfiles
    path_to_candidate_mutation_table = args.candidate_mutation_table
    path_to_cov_mat_raw = args.cov_mat_raw
    path_to_cov_mat_norm = args.cov_mat_norm
    dim=args.dim
    if ('2-case/temp' in path_to_cov_mat_raw) or (path_to_cov_mat_raw == ''):
        # target output is a dummy file and will not be calculated
        Path(path_to_cov_mat_raw).touch()
        flag_cov_raw = False
    else:
        flag_cov_raw = True
    if ('2-case/temp' in path_to_cov_mat_norm) or (path_to_cov_mat_raw == ''):
        # target output is a dummy file and will not be calculated
        Path(path_to_cov_mat_norm).touch()
        flag_cov_norm = False
    else:
        flag_cov_norm = True
    main(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files,
         path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_to_cov_mat_raw, path_to_cov_mat_norm,
         flag_cov_raw, flag_cov_norm,dim)
