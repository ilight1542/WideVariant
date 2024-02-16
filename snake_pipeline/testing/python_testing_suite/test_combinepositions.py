#options for testing:
"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

import sys
import unittest
import os
import glob
import pandas as pd

sys.path.append('./scripts/')
from combine_positions import combine_positions
from combine_positions import generate_positions_snakemake



class TestMyFunction(unittest.TestCase):

    ## set variables to outpaths
    def parse_in_out_variables(self):
        cwd=os.getcwd()
        path_to_samples_csv = glob.glob(f'{cwd}/testing/test_data/sample_csvs/*_samples.csv')
        experiment_names=[path.split('/')[-1].replace('_samples.csv', '') for path in path_to_samples_csv]
        outdirs=[f'testing/test_data/results/{experiment}' for experiment in experiment_names]
        arguments_for_testing=[]
        expected_allpositions_file_list=[]
        for outdir,samples_csv in zip(outdirs,path_to_samples_csv):
            parsed_samples_csv=pd.read_csv(samples_csv, sep=',',header=0)
            cladeIDs=parsed_samples_csv['Group'].unique()
            expected_positions_files_list_of_lists=[]
            for cladeID in cladeIDs:
                expected_allpositions_file_list.append(f'{outdir}/2-case/temp/group_{cladeID}_allpositions.npz')
                arguments_for_testing.append([samples_csv,mapping_file_string,ref_genome,outdir,out_file_string])
                sampleIDs=parsed_samples_csv['Sample'][parsed_samples_csv['Reference'] == ref_genome]
                ref_genomes=parsed_samples_csv['Reference'][parsed_samples_csv['Reference'] == ref_genome]
                outgroup_bools=parsed_samples_csv['Outgroup'][parsed_samples_csv['Reference'] == ref_genome]
                expected_positions_file_list =[]
                for sampleID,ref_genome,outgroup_bool in zip(sampleIDs,ref_genomes,outgroup_bools):
                    expected_positions_file = f'{outdir}/2-case/temp/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.npz'
                    expected_positions_file_list.append(expected_positions_file)
                expected_positions_files_list_of_lists.append(expected_positions_file_list)
            arguments_for_testing.append([expected_positions_file_list,])

    def execute_tests(self,arguments_for_testing,expected_allpositions_file_list,expected_num_samples):
        [positions_files_list] = arguments_for_testing

        combined_positions_obs = generate_positions_snakemake(positions_files_list, REFGENOMEDIRECTORY)
        
        ## check that variants are generated
            # depends on input
        ## check that universal variants are removed
            # depends on input

        with positions_files_list:
            combine_positions(path_to_positions_files, path_to_output_p_file, path_to_outgroup_boolean_file, REFGENOMEDIRECTORY)

        # return self.verificationErrors

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_list,plot_files_list,samples_size_list=self.parse_in_out_variables()
        for args_for_test, expected_plot_files_for_test, expected_num_samples_for_test in zip(args_list,plot_files_list,samples_size_list):
            self.execute_tests(args_for_test, expected_plot_files_for_test, expected_num_samples_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
