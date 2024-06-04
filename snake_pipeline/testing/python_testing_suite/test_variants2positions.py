#options for testing:
"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

## NOTE:
## maxFQ is hardcoded in script and should be changed to a config parameter of snakemake!!

import sys
import unittest
from unittest.mock import patch
import os
import glob
import pandas as pd
import numpy as np
from Bio.Data import IUPACData
from Bio import SeqIO


sys.path.append('../scripts/')
import gus_helper_functions as ghf
from variants2positions import get_fq_score_of_simple_call, generate_positions_single_sample

maxFQ = -30


class TestMyFunction(unittest.TestCase):

    def get_var_pos_raw_dict(self,file, experiment_name):
            var_pos_raw_dict_of_exp = {}
            with open(file, 'r') as fid: 
                header = fid.readline().strip().split()[1:] ## read header without sample column
                for line in fid:
                    line = line.strip().split(',')
                    smplname = line[0]
                    variants_in_sample = line[1:]
                    var_pos_raw_dict_of_exp[f'{experiment_name}_{smplname}'] = {location: int(variant) for location, variant in zip(header, variants_in_sample)}
            return var_pos_raw_dict_of_exp
            
    def parse_in_out_variables(self, maxFQ):
        cwd=os.getcwd()
        results_subdir = 'test_data/results'
        testcases = glob.glob(f'{cwd}/{results_subdir}/*') #/0-used_input_data/samples.csv')
        arguments_for_testing = {}
        expected_output_dict = {}
        var_pos_raw_dicts = {}
        for testcase in testcases:
            experiment_name = testcase.split('/')[-1]
            
            ## read in variants_raw file used for simulation of SNVs to dictionary
            ## get positions which are variable in genome
            var_pos_raw_dicts[experiment_name] = self.get_var_pos_raw_dict(f'{testcase}/0-used_input_data/variant_file.csv', experiment_name)
        
            arguments_for_testing[experiment_name] = []
            expected_output_dict[experiment_name] = []
            samplecsv = pd.read_csv(f'{testcase}/0-used_input_data/samples.csv')
            for sid, sample_row in samplecsv.iterrows():
                sampleID = sample_row["Sample"]
                reference = sample_row["Reference"]
                outgroup = sample_row["Outgroup"] 

                path_to_variant_vcf = f'{testcase}/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz'
                path_to_output_positions = f'{results_subdir}/{experiment_name}/2-case/temp/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.npz'
                path_to_ref = f'{cwd}/{reference}/'
                outgroup_bool = outgroup
                
                arguments_for_testing[experiment_name].append([path_to_variant_vcf,path_to_output_positions,maxFQ,path_to_ref,outgroup_bool])
                ## get genome length, num chromosomes and length of each chromosome
                [chr_starts,genome_length,scafnames] = ghf.genomestats(path_to_ref)
                expected_output_dict[experiment_name].append([sampleID, reference, outgroup, chr_starts, genome_length, scafnames])

        return arguments_for_testing, expected_output_dict, var_pos_raw_dicts

    
    def execute_get_fq_score_of_simple_call_tests(self): 
        ## Test FQ string parsing
        ## generate different input data
        unambiguous_nts = ['A', 'T', 'C', 'G']
        fq_score = -154.5
        vcf_info = f'DP=54;VDB=0.36;SGB=-0.69;MQSB=0.21;MQ0F=0;AF1=1;AC1=2;DP4=0,0,31,9;MQ=41;FQ={fq_score};AC=2;AN=2'
        complex_call = ['A', 'C,T', f'DP=25;VDB=0.823;SGB=-0.692;MQSB=0.988;MQ0F=0;AF1=1;AC1=2;DP4=0,0,18,3;MQ=40;FQ={fq_score};AC=2,0;AN=2'] ## complex call
        missing_FQ_score = ['A', 'C', 'DP=54;VDB=0.36;SGB=-0.69;MQSB=0.21;MQ0F=0;AF1=1;AC1=2;DP4=0,0,31,9;MQ=41;AC=2;AN=2'] ## missing FQ score
        ## for simple calls: 
        for ref_nt in IUPACData.ambiguous_dna_values.keys():
            if ref_nt not in unambiguous_nts:
                ref_is_ambiguous = True
            else:
                ref_is_ambiguous = False
            for alt_nt in IUPACData.ambiguous_dna_values.keys():
                if alt_nt not in unambiguous_nts:
                    alt_is_ambiguous = True
                else:
                    alt_is_ambiguous = False
                simple_call_remove_ambig_fq_score_obs = get_fq_score_of_simple_call(ref_nt, alt_nt, vcf_info, remove_ambigous_call = True)
                if (ref_is_ambiguous | alt_is_ambiguous):
                    ## ref or alt is ambiguous --> no fq score should be reported 
                    with self.subTest(msg=f'test_fq-score_parsing__simple_call_ref{ref_nt}_alt{alt_nt}_remove_ambig'):
                        self.assertTrue( np.isnan(simple_call_remove_ambig_fq_score_obs) )
                else:
                    ## neither ref nor alt is ambiguous --> fq score should be reported 
                    with self.subTest(msg=f'test_fq-score_parsing__simple_call_ref{ref_nt}_alt{alt_nt}_remove_ambig'):
                        self.assertEqual(simple_call_remove_ambig_fq_score_obs, fq_score)
                ## ambiguous sites are not removed --> fq score should be reported 
                simple_call_keep_ambig_fq_score_obs = get_fq_score_of_simple_call(ref_nt, alt_nt, vcf_info, remove_ambigous_call = False)
                with self.subTest(msg=f'test_fq-score_parsing__simple_call_ref{ref_nt}_alt{alt_nt}_keep_ambig'):
                    self.assertEqual(simple_call_keep_ambig_fq_score_obs, fq_score)
        ## complex call --> no fq score should be reported 
        complex_call_fq_score_obs = get_fq_score_of_simple_call(*complex_call)
        with self.subTest(msg=f'test_fq-score_parsing__complex_call'):
            self.assertTrue( np.isnan(complex_call_fq_score_obs) )
        ## Fq score missing --> script should fail (need to be patched for sys.exit(0)!)
        with self.subTest(msg=f'test_fq-score_parsing__missing_FQ_score'):
            with patch('sys.exit') as mock_exit:
                get_fq_score_of_simple_call(*missing_FQ_score)
                mock_exit.assert_called_with(0) # Check if sys.exit was called with the expected return value

    def execute_generate_positions_single_sample_tests(self,arguments_for_testing,expected_output_dict, var_pos_raw_dicts, test_case):
        ## variants2positions
        for single_arg_test, expected_sampleinfo in zip(arguments_for_testing, expected_output_dict):
            path_to_output_positions = single_arg_test[1]
            [sampleID, reference, outgroup, chr_starts, genome_length, scafnames] = expected_sampleinfo
            generate_positions_single_sample(*single_arg_test)
            

            ## check if output_positions file exists        
            with self.subTest(msg=f'{sampleID}_{reference}__expected_positions_exists'):
                self.assertTrue(os.path.exists(path_to_output_positions))

            ## determine size of positions variable
            position_npzfile = np.load(path_to_output_positions)
            npz_file_attribute = position_npzfile.files
            ## check if 'Positions' is present as attribute in input file
            with self.subTest(msg=f'variants2positions_{test_case}__expected_npz_file_attribute'):
                self.assertTrue('Positions' in npz_file_attribute)
            positions_obs = position_npzfile['Positions']
            ## check if not more posotions are in observed positions as reference has nt
            with self.subTest(msg=f'variants2positions_{sampleID}_{reference}__expected_positions_length'):
                self.assertTrue(positions_obs.shape[0]<=genome_length)
            ## check if not more chromosomes are stated in observed positions as reference has chromosomes
            with self.subTest(msg=f'variants2positions_{sampleID}_{reference}__expected_chromosome_number'):
                self.assertTrue(np.unique(positions_obs[:,0]).shape[0] <= len(scafnames))
            ## check for each chromosome that position does not exceed chromosome length
            for chr_id, chr_name in enumerate(scafnames):
                positions_on_chr_idx = (positions_obs[:, 0] == chr_id+1)
                if len(positions_obs[ positions_on_chr_idx, 1 ]) == 0:
                    continue
                positions_on_chromosome = positions_obs[ positions_on_chr_idx, 1 ]
                if chr_id+1 < len(scafnames):
                    chr_length = chr_starts[chr_id+1] - chr_starts[chr_id]
                else:
                    chr_length = genome_length - chr_starts[chr_id]
                with self.subTest(msg=f'variants2positions_{sampleID}_{reference}_{chr_name}__pos_not_exceeding_chromosome_boundaries'):
                    self.assertTrue(max(positions_on_chromosome) <= chr_length)
            
            ## check if outgroup samples have empty position variable
            if outgroup == 1: 
                with self.subTest(msg=f'variants2positions_{sampleID}_{reference}__expected_no_var_position_in_outgroup'):
                    self.assertTrue(positions_obs.shape[0] == 0)
            
            ## check if number of positions and idx in matrix matches expectation
            for location, variant in var_pos_raw_dicts[sampleID].items():
                if variant == 1:
                    chr_var, pos = location.split('_')
                    chr_id = np.where(chr_var == scafnames)[0][0]+1 ## contigs are 1-based!
                    position_exp = np.array([[chr_id, int(pos)]])
                    with self.subTest(msg=f'{test_case}_{sampleID}_{reference}_{location}__expected_variable_pos'):
                        self.assertTrue((position_exp == positions_obs).all(axis = 1).any())

        
        
    #@patch('variants2positions.get_fq_score_of_simple_call')  # Patch sys.exit in function
    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_dict,expected_output_dict, variant_raw_dict = self.parse_in_out_variables(maxFQ)

        self.execute_get_fq_score_of_simple_call_tests()
        for test_case in args_dict.keys():
            self.execute_generate_positions_single_sample_tests(args_dict[test_case], expected_output_dict[test_case], variant_raw_dict[test_case], test_case)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
