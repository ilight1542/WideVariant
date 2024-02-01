#options for testing:
"""
Additional testing script to check if vcf2quals_snakemake.py runs on all test cases properly
"""

import sys
import unittest
import os
import glob
import pandas as pd
from Bio import SeqIO
import numpy as np

sys.path.append('../scripts/')
from vcf2quals_snakemake import round_half_up, vcf_to_quals_snakemake

class TestMyFunction(unittest.TestCase):

    def execute_round_half_up_test(self):
        ## test round_half_up function
        # NOTE: 
        # The function is used because python round() rounds every other .5 val down (to keep the mean equal) 
        for i in range(100):
            with self.subTest(msg=f'Round_up_on_{str(i)}.5'):
                self.assertEqual(round_half_up(i-0.5), i) 
        print('round up half checked\n\n\n')

    ## set variables to outpaths
    def parse_in_out_variables(self):
        cwd=os.getcwd()
        testcases = glob.glob(f'{cwd}/test_data/results/*') #/0-used_input_data/samples.csv')
        arguments_for_testing = {}
        var_raw_dicts = {}
        for testcase in testcases:
            experiment_name = testcase.split('/')[-1]
            ## read in variants_raw file to dictionayr
            var_raw_dicts[experiment_name] = {}
            arguments_for_testing[experiment_name] = []
            with open(f'{testcase}/0-used_input_data/variant_file.csv', 'r') as fid: 
                for lineid, line in enumerate(fid):
                    line = line.strip().split(',')
                    if lineid == 0:
                        header = line
                    else:
                        var_raw_dicts[experiment_name][line[0]] = {location: int(variant) for location, variant in zip(header[1:], line[1:])}
            samplecsv = pd.read_csv(f'{testcase}/0-used_input_data/samples.csv')
            for sid, sample_row in samplecsv.iterrows():
                sampleid = sample_row["Sample"]
                reference = sample_row["Reference"]
                outgroup = sample_row["Outgroup"]
                path_to_vcf_files = f'{testcase}/1-mapping/vcf/{sampleid}_ref_{reference}_aligned.sorted.strain.vcf.gz'
                path_to_ref = f'{cwd}/{reference}/'
                path_tp_outfiles = f'test_data/results/{experiment_name}/1-mapping/vcf/{sampleid}_ref_{reference}_outgroup{outgroup}.quals.npz' ## need to be relative path!
                arguments_for_testing[experiment_name].append([path_to_vcf_files, path_tp_outfiles, path_to_ref])
        return arguments_for_testing, var_raw_dicts
            
    def execute_vcf2quals_test(self,arguments_for_testing,var_raw_dict,test_case):
        
        for single_arg_test in arguments_for_testing:
            path_to_vcf_files, path_to_quals_file, reference_path = single_arg_test
            vcf_to_quals_snakemake(*single_arg_test)
        
            reference_name = reference_path.split('/')[-1]
            ## get length of genome
            contig_length_dict = {}
            with open(f'{reference_path}/genome.fasta', 'r') as fasta:
                for record in SeqIO.parse(fasta, 'fasta'):
                    contig_length_dict[record.id] = len(record)
            genome_length = sum(contig_length_dict.values())
            ## get sample name
            smpl = path_to_quals_file.split('/')[-1].split('_ref')[0]
            ## read in quals file
            quals = np.load(path_to_quals_file)['quals']
            ## check if quals file is length of reference
            quals_length = np.shape(quals)[0]
            with self.subTest(msg=f'{smpl}_{reference_name}__expected_num_bases'):
                self.assertEqual(quals_length, genome_length) ## note for each sample we read 2 lines!
            ## check where smpl has variant and whether variant has qual score 
            for location, variant in var_raw_dict.items():
                if variant == 1:
                    contig_var, pos = location.split('_')
                    pos_in_genome = 0
                    for contig_genome, contig_length in genome_length.items():
                        if contig_genome != contig_var:
                            pos_in_genome+=contig_length
                        else:
                            break
                    var_pos_genome = pos_in_genome+int(pos)-1 ## get variant position within genome
                    with self.subTest(msg=f'{smpl}_{reference_name}_{location}__expected_qual_score'):
                        self.assertTrue((quals[var_pos_genome] == [0])[0])
        print(f'vcf2quals checked on {test_case}\n\n\n')

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        self.execute_round_half_up_test()
        args_dict,variant_raw_dict = self.parse_in_out_variables()

        for test_case in args_dict.keys():
            self.execute_vcf2quals_test(args_dict[test_case], variant_raw_dict[test_case], test_case)
        

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31

