#options for testing:
"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

import sys
import unittest
import os
import glob
import pandas as pd

sys.path.append(f'../scripts/')
from bowtie2qc import plot_bowtieqc



class TestMyFunction(unittest.TestCase):

    ## set variables to outpaths
    def parse_in_out_variables(self):
        cwd=os.getcwd()
        testcase_dirs = glob.glob(f'{cwd}/test_data/results/*')
        arguments_for_testing_dict = {}
        expected_plot_files_dict = {}
        expected_num_samples_dict = {}
        expected_text_file_dict = {}
        for testcase_dir in testcase_dirs:
            experiment_name = testcase_dir.split('/')[-1]
            arguments_for_testing_dict[experiment_name] = []
            expected_plot_files_dict[experiment_name] = []
            expected_num_samples_dict[experiment_name] = []
            expected_text_file_dict[experiment_name] = []
            outdir=f'{testcase_dir}/1-mapping/bowtie2_qc'
            samples_csv_path=f'{testcase_dir}/0-used_input_data/samples.csv'
            samples_csv=pd.read_csv(f'{testcase_dir}/0-used_input_data/samples.csv', sep=',',header=0)
            ref_genomes=samples_csv['Reference'].unique()
            for ref_genome in ref_genomes:
                expected_num_samples = len(samples_csv[samples_csv['Reference'] == ref_genome])
                expected_num_samples_dict[experiment_name].append(expected_num_samples)
                expected_plot_files=[f'{outdir}/alignment_stats_ref_{ref_genome}_number_aligned_once_histogram.png',
                                    f'{outdir}/alignment_stats_ref_{ref_genome}_percent_aligned_once_histogram.png',
                                    f'{outdir}/alignment_stats_ref_{ref_genome}_overall_alignment_histogram.png']
                mapping_file_string=f'{testcase_dir}/1-mapping/bowtie2/'
                out_file_string=f'alignment_stats_ref_{ref_genome}'
                expected_text_file=f'{outdir}/{out_file_string}.txt'
                arguments_for_testing_dict[experiment_name].append([samples_csv_path,mapping_file_string,ref_genome,outdir,out_file_string])
                expected_plot_files_dict[experiment_name].append(expected_plot_files)
                expected_text_file_dict[experiment_name].append(expected_text_file)
        return arguments_for_testing_dict,expected_text_file_dict,expected_plot_files_dict,expected_num_samples_dict

    def execute_tests(self,arguments_for_testing,expected_text_file,expected_plot_files,expected_num_samples):
        [_,_,ref_genome,outdir,_] = arguments_for_testing
        experiment_name = outdir.split('/')[-3]
        ## run plotting to compare
        plot_bowtieqc(*arguments_for_testing)
        ## check if mapping stats grep file exists 
        with self.subTest(msg=f'{experiment_name}_{ref_genome}__text_file'):
            self.assertTrue(os.path.exists(expected_text_file))
        ## check number of lines after grep argument of QC script
        with open(expected_text_file, 'r') as fid:
            num_lines_obs = len(fid.readlines())
        ## assess if all samples are in collapsed mapping style file
        with self.subTest(msg=f'{experiment_name}_{ref_genome}__expected_num_samples'):
            self.assertEqual(num_lines_obs, expected_num_samples * 2) ## note for each sample we read 2 lines!
        ## assess if all plots are generated
        for i, expected_file in enumerate(expected_plot_files):
            with self.subTest(msg=f'{experiment_name}_{ref_genome}__expected_plot{i}'):
                self.assertTrue(os.path.exists(expected_file))
        # return self.verificationErrors

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_dict,text_file_dict,plot_files_dict,samples_size_dict=self.parse_in_out_variables()
        for testcase in args_dict.keys():
            for args_for_test, expected_text_files_for_test, expected_plot_files_for_test, expected_num_samples_for_test in zip(args_dict[testcase],text_file_dict[testcase],plot_files_dict[testcase],samples_size_dict[testcase]):
                self.execute_tests(args_for_test, expected_text_files_for_test, expected_plot_files_for_test, expected_num_samples_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31

