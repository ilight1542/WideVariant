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
import gus_helper_functions as gus

class TestMyFunction(unittest.TestCase):
    def test_read_fasta(self):
        self.assertRaises(gus.read_fasta('testing/test_data/gus_test_data/non_findable_genome'))
        self.assertIsInstance(read_in_fasta, SeqIO.Seq)

    def test_genomestats(self):
        pass

    def test_p2chrpos(self):
        pass


    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_list,plot_files_list,samples_size_list=self.parse_in_out_variables()
        for args_for_test, expected_plot_files_for_test, expected_num_samples_for_test in zip(args_list,plot_files_list,samples_size_list):
            [_,_,ref_genome,outdir,out_file_string] = args_for_test
            experiment_name = outdir.split('/')[-1]
            with self.subTest(msg=f'TESTING: {experiment_name}_{ref_genome}'):
                self.execute_tests(args_for_test, expected_plot_files_for_test, expected_num_samples_for_test)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
