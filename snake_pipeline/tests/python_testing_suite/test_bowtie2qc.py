#options for testing:

"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

import unittest
import os
from ../scripts/etc/bowtie2qc import main

## set variables to outpaths
path_to_samples = '../../samples.csv'
path_to_mappingsstats = 'results/1-mapping/bowtie2'
reference_genome = ''
outdir = '../../results'
out_file_string = f'results/1-mapping/bowtie2_qc/alignment_stats_ref_{reference_genome}'

## set arguments for testing
arguments_no_outgroup = ['path_to_samples_no_outgroup', 'path_to_mappingsstats_no_outgroup', reference_genome, outdir, out_file_string]
arguments_all_outgroup = ['path_to_samples_all_outgroup', 'path_to_mappingsstats_no_outgroup', reference_genome, outdir, out_file_string]
arguments_missing_smpl = ['path_to_samples_missing_smpl', 'path_to_mappingsstats_missing_smpl', reference_genome, outdir, out_file_string]
expected_files = [f'{outdir}/{out_file_string}_number_aligned_once_histogram.png',
                  f'{outdir}/{out_file_string}_percent_aligned_once_histogram.png',
                  f'{outdir}/{out_file_string}_overall_alignment_histogram.png'] 

class TestMyFunction(unittest.TestCase):
    def test_outputs(self):
        ## Read in number of samples expected
        num_samples = 0
        with open('../../samples.csv', 'r') as fid:
            for line in fid:
                if reference_genome in line:
                    num_samples += 1
        ## check if mapping stats grep file exists 
        self.assertTrue(os.path.exists(f'{outdir}/{out_file_string}.txt'))
        ## check number of lines after grep argument of QC script
        with open(f'{outdir}/{out_file_string}.txt', 'r'):
            num_samples_obs = len(fid.readlines())
        ## assess if all samples are in collapsed mapping style file
        self.assertEqual(num_samples_obs, num_samples * 2) ## note for each sample we read 2 lines!
        ## assess if all plots are generated
        for expected_file in expected_files:
            self.assertTrue(os.path.exists(expected_file))
    
    ## Check outputs for run without outgroup
    def test_no_outgroup(self):
        ## run script on samples
        main(arguments_no_outgroup)
        self.test_outputs()

    ## Check outputs when all samples are outgroup
    def test_all_outgroup(self):
        ## run script on samples
        main(arguments_all_outgroup)
        self.test_outputs()

    ## Check outputs for run when some files are missing after mapping
    def test_missing_bowtie_file_for_one_sample(self):
        ## run script on samples
        main(arguments_missing_smpl)
        self.test_outputs()

if __name__ == '__main__':
    unittest.main()