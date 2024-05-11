
"""
Additional testing script to check if combine_positions.py runs on all test cases properly
"""

import sys
import unittest
import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.append('../scripts/')
import gus_helper_functions as ghf
from combine_positions import combine_positions, generate_positions_snakemake


class TestMyFunction(unittest.TestCase):

    def get_var_pos_ingroup_raw_dict(self,file, ingroup_samples, header = True):
        var_pos_ingroup_raw_of_exp_dict = {}
        with open(file, 'r') as fid: 
            header = fid.readline().strip().split(',')[1:] ## read header 
            data = [line.strip().split(',') for line in fid.readlines()] ## read remaining data 
            ingroup_data = np.array([line[1:] for line in data if line[0] in ingroup_samples]) ## read rest of data without sample names if sample an ingroup sample
            ingroup_data = ingroup_data.astype(int)
        if ingroup_data.size > 0:
            num_samples = np.shape(ingroup_data)[0]
            for pos_idx in range(np.shape(ingroup_data)[1]):
                num_variable_pos = len(np.where(ingroup_data[:,pos_idx] == 1)[0])
                var_pos_ingroup_raw_of_exp_dict[header[pos_idx]] = ((num_variable_pos>0) & (num_variable_pos<num_samples)) ## if more than 1 entry in data column & less then all samples have entry --> variable position
        else:
            var_pos_ingroup_raw_of_exp_dict = {location: False for location in header} ## all given positions are invariable or no ingroup sample is given
        return var_pos_ingroup_raw_of_exp_dict

    def update_moved_files_in_snakemake_generated_txtfiles(self, infile, outfile, old_path_string, new_path_string):
        with open(infile, 'r') as old_pathnames_file, open(outfile, 'w') as updated_pathnames_out:
            for line in old_pathnames_file:
                # Replace the old path with the new path
                modified_line = line.replace(old_path_string, new_path_string)
                updated_pathnames_out.write(modified_line)

    def parse_in_out_variables(self):
        cwd=os.getcwd()
        results_subdir = 'test_data/results'
        testcases = glob.glob(f'{cwd}/{results_subdir}/*') #/0-used_input_data/samples.csv')
        arguments_for_testing = {}
        var_pos_ingroup_raw_dicts = {}
        smpl_list = {}
        for testcase in testcases:
            experiment_name = testcase.split('/')[-1]
            ## read in variants_raw file used for simulation of SNVs to dictionary
            samplecsv = pd.read_csv(f'{testcase}/0-used_input_data/samples.csv')
            ingroup_samples = samplecsv.loc[samplecsv['Outgroup'] == 0, 'Sample'].str.replace(f'{experiment_name}_', '').values

            ## get positions which are variable in genome
            var_pos_ingroup_raw_dicts[experiment_name] = self.get_var_pos_ingroup_raw_dict(f'{testcase}/0-used_input_data/variant_file.csv', ingroup_samples)
            ## get unique sets of reference and group
            samplecsv_uniq_ref_grp_pairs = samplecsv[['Reference', 'Group']].drop_duplicates()
            arguments_for_testing[experiment_name] = []
            smpl_list[experiment_name] = []
            for sid, sample_row in samplecsv_uniq_ref_grp_pairs.iterrows():
                reference = sample_row["Reference"]
                groupid = sample_row["Group"] 
                
                path_to_called_pos_file = f'{testcase}/2-case/temp/group_{groupid}_string_file_other_p_to_consider.txt'
                path_to_called_pos_moved_paths_file = f'{testcase}/2-case/temp/group_{groupid}_string_file_other_p_to_consider_moved_filenames.txt'
                ## modify in file the path names as they have been moved!
                self.update_moved_files_in_snakemake_generated_txtfiles(path_to_called_pos_file, path_to_called_pos_moved_paths_file, 
                                                                        old_path_string = 'results/2-case',
                                                                        new_path_string = f'test_data/results/{experiment_name}/2-case')
                
                path_to_outgroup_file = f'{testcase}/2-case/temp/group_{groupid}_string_outgroup_bool.txt'
                path_to_ref = f'{cwd}/{reference}/'
                path_to_outfile = f'{results_subdir}/{experiment_name}/2-case/temp/group_{groupid}_allpositions.npz' ## need to be relative path!
                arguments_for_testing[experiment_name].append([path_to_called_pos_moved_paths_file, path_to_outfile, path_to_outgroup_file, path_to_ref])
                
                smpls_in_groupid = samplecsv.loc[(samplecsv['Reference'] == reference) & (samplecsv['Group'] == groupid), 'Sample'].values ## contains per group id entry the list of samples associated with it --> needed to check later! 
                smpl_list[experiment_name].append(smpls_in_groupid) 

        return arguments_for_testing, var_pos_ingroup_raw_dicts, smpl_list

    def get_ingroup_positions_files(self, path_to_called_pos_moved_paths_file, path_to_outgroup_file):
        ## get just called position files which are in ingroup
        positions_files_lists = []
        with open(path_to_called_pos_moved_paths_file) as called_pos_file, open(path_to_outgroup_file) as outgroup_file:
            for called_pos_line, outgroup_line in zip(called_pos_file, outgroup_file):
                if int(outgroup_line.strip()) == 0:
                    positions_files_lists.append(called_pos_line.strip())
        return positions_files_lists

    def execute_combine_positions_prep_tests(self,arguments_for_testing, smpl_list, test_case):
        ## This tests the snakemake rule if all files are taken as input for combine positions
        for single_arg_test, smpl_list_of_grpid in zip(arguments_for_testing, smpl_list):
            path_to_called_pos_moved_paths_file, _, _, _ = single_arg_test
            smplids_combine_positions_prep_file = []
            with open(path_to_called_pos_moved_paths_file, 'r') as fid: 
                for line in fid:
                    line = line.strip().split('/')[-1].split('_ref_')[0]
                    smplids_combine_positions_prep_file.append(line.split('_ref_')[0])
            ## check if sample ids expected for the group and the sample ids given in file for that group are the same length 
            with self.subTest(msg=f'combine_positions_prep_{test_case}__expected_num_sampleids'):
                self.assertEqual(len(smplids_combine_positions_prep_file), len(smpl_list_of_grpid))
            ## check further if both have the same sampleids 
            for exp_smplid in smpl_list_of_grpid:
                with self.subTest(msg=f'combine_positions_prep_{test_case}__expected_sampleids_used'):
                    self.assertTrue(exp_smplid in smplids_combine_positions_prep_file)
        print(f'combine_positions_prep checked on {test_case}\n\n\n')          

    def execute_generate_positions_snakemake_tests(self,arguments_for_testing, var_pos_ingroup_raw_dicts, test_case):
        
        ## check if generate_positions_snakemake runs fine if path_to_called_pos_moved_paths_file is empty and returns empty string
        path_to_ref = arguments_for_testing[0][3] ## get first path to ref in list 
        combined_pos_empty_input = generate_positions_snakemake([], path_to_ref)
        with self.subTest(msg=f'hardcoded_test__empty_input_lists_generate_positions_snakemake'):
            self.assertTrue(len(combined_pos_empty_input) == 0)

        for single_arg_test in arguments_for_testing:
            path_to_called_pos_moved_paths_file, _, path_to_outgroup_file, path_to_ref = single_arg_test
            ## get files which are in ingroup only
            positions_files_list = self.get_ingroup_positions_files(path_to_called_pos_moved_paths_file, path_to_outgroup_file)
            ## check if files in path_to_called_pos_moved_paths_file are present and contain correct 'Positions' attribute
            for fileid, file in enumerate(positions_files_list):
                ## check if file is present
                with self.subTest(msg=f'generate_positions_snakemake_{test_case}__expected_positions_file{fileid}'):
                    self.assertTrue(os.path.exists(file))
                npz_file_attribute = np.load(file, mmap_mode='r').files
                ## check if Positions is present in input file
                with self.subTest(msg=f'generate_positions_snakemake_{test_case}__expected_positions_file_attribute{fileid}'):
                    self.assertTrue('Positions' in npz_file_attribute)
            combined_pos = generate_positions_snakemake(positions_files_list, path_to_ref)
            ## get length of genome
            [chr_starts,genome_length,scafnames] = ghf.genomestats(path_to_ref)
            ## get set variable variants 
            for location, is_variable_pos in var_pos_ingroup_raw_dicts.items():
                contig_var, pos = location.split('_')
                contig_id = np.where(scafnames == contig_var)[0][0] + 1
                var_pos_genome = ghf.chrpos2index(np.array( ([contig_id, int(pos)], ) ), chr_starts) ## get variant position within genome
                if is_variable_pos:
                    ## check if variable position in ingroup samples present in combined position
                    with self.subTest(msg=f'generate_positions_snakemake_{test_case}_{location}__expected_variable_pos_ingroup'):
                        self.assertTrue(var_pos_genome in combined_pos)
                else:
                    ## check if invariable position in ingroup samples not in combined position
                    with self.subTest(msg=f'generate_positions_snakemake_{test_case}_{location}__unexpected_invariable_pos_ingroup'):
                        self.assertTrue(var_pos_genome not in combined_pos)
        print(f'generate_positions_snakemake checked on {test_case}\n\n\n')        

    def execute_combine_positions_tests(self,arguments_for_testing, var_pos_ingroup_raw_dicts, test_case):
        
        for single_arg_test in arguments_for_testing:
            combine_positions(*single_arg_test)
            path_to_called_pos_moved_paths_file, path_to_outfile, path_to_outgroup_file, path_to_ref = single_arg_test
            
            ## check if npz file is present
            with self.subTest(msg=f'combine_positions_{test_case}__expected_npz_file'):
                self.assertTrue(os.path.exists(path_to_outfile))
            
            combined_positions_outfile = np.load(path_to_outfile, mmap_mode='r')
            npz_file_attribute = combined_positions_outfile.files
            ## check if 'p' is present as attribute in input file
            with self.subTest(msg=f'combine_positions_{test_case}__expected_npz_file_attribute'):
                self.assertTrue('p' in npz_file_attribute)
            
            ## load positions from npz file
            combined_positions = combined_positions_outfile['p']
            
            ## get starts of genome
            [chr_starts,_,scafnames] = ghf.genomestats(path_to_ref)
            ## get set variable variants 
            for location, is_variable_pos in var_pos_ingroup_raw_dicts.items():
                contig_var, pos = location.split('_')
                contig_id = np.where(scafnames == contig_var)[0][0] + 1
                var_pos_genome = ghf.chrpos2index(np.array( ([contig_id, int(pos)], ) ), chr_starts) ## get variant position within genome
                if is_variable_pos:
                    ## check if variable position in ingroup samples present in combined position
                    with self.subTest(msg=f'combine_positions_{test_case}_{location}__expected_variable_pos_ingroup'):
                        self.assertTrue(var_pos_genome in combined_positions)
                else:
                    ## check if invariable position in ingroup samples not in combined position
                    with self.subTest(msg=f'combine_positions_{test_case}_{location}__unexpected_invariable_pos_ingroup'):
                        self.assertTrue(var_pos_genome not in combined_positions)
        print(f'combine_positions checked on {test_case}\n\n\n')       

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_dict,variant_raw_dict, smpl_list = self.parse_in_out_variables()

        for test_case in args_dict.keys():
            ## test each function of combine positions script independently 
            self.execute_combine_positions_prep_tests(args_dict[test_case], smpl_list[test_case], test_case) ## NOTE: This script tests if the snakemake input was correctly run before
            self.execute_generate_positions_snakemake_tests(args_dict[test_case], variant_raw_dict[test_case], test_case)
            self.execute_combine_positions_tests(args_dict[test_case], variant_raw_dict[test_case], test_case)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
