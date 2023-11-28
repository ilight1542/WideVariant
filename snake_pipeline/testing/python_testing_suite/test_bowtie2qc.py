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
from bowtie2qc import plot_bowtieqc



class TestMyFunction(unittest.TestCase):

    ## set variables to outpaths
    def parse_in_out_variables(self):
        cwd=os.getcwd()
        path_to_samples_csv = glob.glob(f'{cwd}/testing/test_data/sample_csvs/*_samples.csv')
        experiment_names=[path.split('/')[-1].replace('_samples.csv', '') for path in path_to_samples_csv]
        outdirs=[f'testing/test_data/results/{experiment}' for experiment in experiment_names]
        arguments_for_testing=[]
        expected_plot_files_list=[]
        expected_num_samples_list=[]
        for outdir,samples_csv in zip(outdirs,path_to_samples_csv):
            parsed_samples_csv=pd.read_csv(samples_csv, sep=',',header=0)
            ref_genomes=parsed_samples_csv['Reference'].unique()
            for ref_genome in ref_genomes:
                expected_num_samples = len(parsed_samples_csv[parsed_samples_csv['Reference'] == ref_genome])
                expected_num_samples_list.append(expected_num_samples)
                expected_plot_files=[f'{outdir}/1-mapping/bowtie2_qc/alignment_stats_ref_{ref_genome}_number_aligned_once_histogram.png',
                                    f'{outdir}/1-mapping/bowtie2_qc/alignment_stats_ref_{ref_genome}_percent_aligned_once_histogram.png',
                                    f'{outdir}/1-mapping/bowtie2_qc/alignment_stats_ref_{ref_genome}_overall_alignment_histogram.png']
                mapping_file_string=f'{outdir}/1-mapping/bowtie2/'
                out_file_string=f'alignment_stats_ref_{ref_genome}'
                arguments_for_testing.append([samples_csv,mapping_file_string,ref_genome,outdir,out_file_string])
                expected_plot_files_list.append(expected_plot_files)
        return arguments_for_testing,expected_plot_files_list,expected_num_samples_list

    def execute_tests(self,arguments_for_testing,expected_plot_files,expected_num_samples):
        [_,_,ref_genome,outdir,out_file_string] = arguments_for_testing
        expected_text_file=f'{outdir}/{out_file_string}.txt'
        experiment_name = outdir.split('/')[-1]
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
        args_list,plot_files_list,samples_size_list=self.parse_in_out_variables()
        for args_for_test, expected_plot_files_for_test, expected_num_samples_for_test in zip(args_list,plot_files_list,samples_size_list):
            self.execute_tests(args_for_test, expected_plot_files_for_test, expected_num_samples_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
