#options for testing:
"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

import sys
import unittest
import os
import glob
import pandas as pd
import numpy as np


sys.path.append('./scripts/')
import gus_helper_functions as ghf
from build_candiate_mutation_table import build_CMT



class TestMyFunction(unittest.TestCase):

    ## set variables to outpaths
    def parse_in_out_variables(self):
        cwd=os.getcwd()
        path_to_samples_csv = glob.glob(f'{cwd}/testing/test_data/sample_csvs/*_samples.csv')
        testing_directory = f'testing/test_data/'
        experiment_names=[path.split('/')[-1].replace('_samples.csv', '') for path in path_to_samples_csv]
        outdirs=[f'testing/test_data/results/{experiment}' for experiment in experiment_names]
        arguments_for_testing=[]
        expected_positions_files_list=[]
        expected_sample_info=[]

        for outdir,samples_csv in zip(outdirs,path_to_samples_csv):
            parsed_samples_csv=pd.read_csv(samples_csv, sep=',',header=0)
            sampleIDs=parsed_samples_csv['Sample']
            ref_genomes=parsed_samples_csv['Reference']
            outgroup_bools=parsed_samples_csv['Outgroup']
            for sampleID,ref_genome,outgroup_bool in zip(sampleIDs,ref_genomes,outgroup_bools):


                build_CMT(path_to_p_file, path_to_sample_names_file, path_to_outgroup_boolean_file, path_to_list_of_quals_files,
                 path_to_list_of_diversity_files, path_to_candidate_mutation_table, path_to_cov_mat_raw, path_to_cov_mat_norm,
                 flag_cov_raw, flag_cov_norm, dim)


                
                

                arguments_for_testing.append([expected_variant_vcf,expected_positions_file,maxFQ,ref_genome_directory,outgroup_bool])
                expected_positions_files_list.append(expected_positions_file)
                expected_sample_info.append([sampleID, ref_genome, outdir, max_expected_length, max_expected_chrom])
        return arguments_for_testing,expected_positions_files_list, expected_sample_info
    def execute_tests(self,arguments_for_testing,expected_positions_file, expected_sample_info):
        ## build_candidate_mutation_table
        build_CMT(*arguments_for_testing)

        ## check if CMT variables are the right shape

        #quals: num_samples x num_pos
        #p: num_pos
        #counts: num_samples x num_pos x 8
        #in_outgroup: num_samples 
        #sampleNames: num_samples
        #indel_counter: num_samples x num_pos x 2

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_list,positions_files_list,sample_info_list=self.parse_in_out_variables()
        for args_for_test, expected_position_files_for_test, expected_sample_info_for_test in zip(args_list,positions_files_list,sample_infos):
            self.execute_tests(args_for_test, expected_position_files_for_test, expected_sample_info_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
