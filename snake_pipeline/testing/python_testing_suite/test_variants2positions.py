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
from variants2positions import generate_positions_single_sample


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
                expected_variant_vcf = f'{outdir}/1-mapping/vcf/{sampleID}_ref_{ref_genome}_aligned.sorted.strain.variant.vcf.gz'
                expected_positions_file = f'{outdir}/2-case/temp/{sampleID}_ref_{ref_genome}_outgroup{outgroup_bool}_positions.npz'
                maxFQ=30
                ref_genome_directory = f'testing/{ref_genome}/'
                [chr_starts,genome_length,scaf_names] = ghf.genomestats(ref_genome_directory)
                max_expected_length = genome_length 
                max_expected_chrom = chr_starts.size
                arguments_for_testing.append([expected_variant_vcf,expected_positions_file,maxFQ,ref_genome_directory,outgroup_bool])
                expected_positions_files_list.append(expected_positions_file)
                expected_sample_info.append([sampleID, ref_genome, outdir, max_expected_length, max_expected_chrom])
        return arguments_for_testing,expected_positions_files_list, expected_sample_info
    def execute_tests(self,arguments_for_testing,expected_positions_file, expected_sample_info):
        ## variants2positions
        [path_to_variant_vcf, path_to_output_positions, maxFQ, ref_genome_directory, outgroup_bool] = arguments_for_testing
        [sampleID, ref_genome, outdir, max_expected_length, max_expected_chrom] = expected_sample_info
        generate_positions_single_sample(*arguments_for_testing)

        ## check if output_positions file exists        
        with self.subTest(msg=f'{sampleID}_{ref_genome}__expected_positions_exists'):
            self.assertTrue(os.path.exists(expected_positions_file))

        ## determine size of positions variable
        position_file = np.load(expected_positions_file)
        positions_obs = position_file['Positions']
        with self.subTest(msg=f'{sampleID}_{ref_genome}__expected_positions_length'):
            self.assertTrue(positions_obs.shape[0]<=max_expected_length)
        with self.subTest(msg=f'{sampleID}_{ref_genome}__expected_chromosome_number'):
            self.assertTrue(max(positions_obs[:,0])<=max_expected_chrom)
        ## check if outgroup samples have empty position variable
        # captured in above test 

        ## check if number of positions matches expectation

        ## based on wgsim with high enough cov/qual/etc

        

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_list,positions_files_list,sample_infos=self.parse_in_out_variables()
        for args_for_test, expected_position_files_for_test, expected_sample_info_for_test in zip(args_list,positions_files_list,sample_infos):
            self.execute_tests(args_for_test, expected_position_files_for_test, expected_sample_info_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
