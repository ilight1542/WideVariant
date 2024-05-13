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
from Bio import SeqIO

sys.path.append('../../scripts/')
import gus_helper_functions as ghf
import pileup2diversity as p2dh

class TestMyFunction(unittest.TestCase):
    def test_generate_nts_ascii(self):
        nt_order=ghf.get_nt_order()
        ord_conversion={'A': 65, 'T': 84, 'C': 67, 'G': 71, 'a': 97, 't': 116, 'c': 99, 'g': 103}
        self.assertEqual(p2dh.generate_nts_ascii()[:4], [ord_conversion[x] for x in ['A','T','C','G']])
        self.assertEqual(p2dh.generate_nts_ascii()[4:], [ord_conversion[x] for x in ['a','t','c','g']])
    
    def test_clean_calls_from_start_and_end(self):
        no_start_or_end=np.fromstring('....,,,,', dtype=np.int8)
        start_fwd_rev_entries = '^X.^X,'
        end_fwd_rev_entries = '$.$,'
        all_start=np.fromstring(start_fwd_rev_entries*4, dtype=np.int8)
        all_end=np.fromstring(end_fwd_rev_entries*4, dtype=np.int8)
        mixed_start_end=np.fromstring(f'{start_fwd_rev_entries}{end_fwd_rev_entries}'*2, dtype=np.int8)
        mixed_start_end_mid=np.fromstring(f'{start_fwd_rev_entries}{end_fwd_rev_entries}..,,'*2, dtype=np.int8)

        all_start_exp_arr = np.tile(np.array([-1, -1, 46, -1, -1, 44], dtype=np.int8), 4)
        all_end_exp_arr = np.tile(np.array([-1, 46, -1, 44], dtype=np.int8), 4)
        mixed_start_end_exp_arr = np.tile(np.array([-1, -1, 46, -1, -1, 44, -1, 46, -1, 44], dtype=np.int8), 2)
        mixed_start_end_mid_exp_arr = np.tile(np.array([-1, -1, 46, -1, -1, 44, -1, 46, -1, 44, 46, 46, 44, 44], dtype=np.int8), 2)

        # test that lengths are correct for clean calls
        self.assertTrue( np.array_equal(p2dh.clean_calls_from_start_and_end(no_start_or_end.copy())[0], no_start_or_end) ) ## check if calls have not been modified
        self.assertEqual( p2dh.clean_calls_from_start_and_end(no_start_or_end.copy())[1], 0) ## no end or start reads have been identified
        self.assertTrue( np.array_equal(p2dh.clean_calls_from_start_and_end(all_start.copy())[0], all_start_exp_arr) ) ## check if all start chars '^' + mapping value ('X') have been flagged
        self.assertEqual( p2dh.clean_calls_from_start_and_end(all_start.copy())[1], 8) ## all are start reads
        self.assertTrue( np.array_equal(p2dh.clean_calls_from_start_and_end(all_end.copy())[0], all_end_exp_arr) ) ## check if all end chars '$' have been flagged
        self.assertEqual( p2dh.clean_calls_from_start_and_end(all_end.copy())[1], 8) ## all are end reads
        self.assertTrue( np.array_equal(p2dh.clean_calls_from_start_and_end(mixed_start_end.copy())[0], mixed_start_end_exp_arr) ) ## check if all start and end chars have been flagged
        self.assertEqual( p2dh.clean_calls_from_start_and_end(mixed_start_end.copy())[1], 8) ## all are either end or start reads
        self.assertTrue( np.array_equal(p2dh.clean_calls_from_start_and_end(mixed_start_end_mid.copy())[0], mixed_start_end_mid_exp_arr) ) ## check if all start and end chars have been flagged and normal reads retained 
        self.assertEqual( p2dh.clean_calls_from_start_and_end(mixed_start_end_mid.copy())[1], 8) ## all are either end or start reads


    def test_parse_indels_into_data(self):
        #TODO finish implementation (checking updating of both calls and data after parse_indels_into_data runs)
        # indel at start of contig
        # indel near start (overlapping based on indelregion to mark +/- #positions as being in indel proximity)
        # indel at end of contig
        ## generate dataframe to update
        test_data=np.zeros((10,40),dtype=int) # 10bp genome length
        
        ## different types of indel entries for calls
        no_start_or_end=np.fromstring('....,,,,', dtype=np.int8)
        all_insertions=np.fromstring('.+1A.+1A.+1A.+1A.+1A.+1A.+1A.+1A', dtype=np.int8)
        two_base_insertion=np.fromstring('.+2AA.+2AA.+2AA.+2AA.+2AA.+2AA.+2AA.+2AA', dtype=np.int8)
        all_deletions=np.fromstring('.-1A.-1A.-1A.-1A.-1A.-1A.-1A.-1A', dtype=np.int8)
        two_base_deletions=np.fromstring('.-2AA.-2AA.-2AA.-2AA.-2AA.-2AA.-2AA.-2AA', dtype=np.int8)

        ## build sets of calls:
        calls_no_indels=np.array([no_start_or_end]*10)
        calls_insertion_start=
        calls_deletion_start=
        calls_twobp_insertion_middle=

        # TODO: check calls get updated as -1 correctly for strings that arent ./,/ATCG

        self.assertEqual(p2dh.parse_indels_into_data(calls_no_indels,test_data),8)
        self.assertEqual(p2dh.parse_indels_into_data(all_end),8)
        self.assertEqual(p2dh.parse_indels_into_data(mixed_start_end),8)

    def test_parse_calls_into_simplecalls(self):
        #TODO check if it works, possibly add more edge cases
        all_ref=np.fromstring('....,,,,', dtype=np.int8)
        no_ref=np.fromstring('TTTTtttt', dtype=np.int8)
        all_ref_removed_bases=np.array(list(np.fromstring('....,,,,', dtype=np.int8))+[-1,-1,-1,-1])
        ref_idx = [0, 1, 2, 3]
        nts_ascii=p2dh.generate_nts_ascii()
        pileup_line = 'pseudo-pileup line'
        pileup_lineid = 0
        
        # test that lengths are correct for clean calls
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(all_ref.copy(), ref_idx[0], nts_ascii, pileup_line, pileup_lineid), 
                                        np.full(8, nts_ascii[0], dtype=np.int8) )) ## ref is 'A' & called ref --> 8xASCII char '65'
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(all_ref.copy(), ref_idx[1], nts_ascii, pileup_line, pileup_lineid), 
                                        np.full(8, nts_ascii[1], dtype=np.int8) )) ## ref is 'T' & called ref --> 8xASCII char '84'
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(all_ref.copy(), ref_idx[2], nts_ascii, pileup_line, pileup_lineid), 
                                        np.full(8, nts_ascii[2], dtype=np.int8) )) ## ref is 'C' & called ref --> 8xASCII char '67'
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(all_ref.copy(), ref_idx[3], nts_ascii, pileup_line, pileup_lineid), 
                                        np.full(8, nts_ascii[3], dtype=np.int8) )) ## ref is 'G' & called ref --> 8xASCII char '71'
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(no_ref.copy(), ref_idx[0], nts_ascii, pileup_line, pileup_lineid),  
                                        np.array([ 84,  84,  84,  84, 116, 116, 116, 116], dtype=np.int8) )) ## ref 'T' and 'A' called  --> 8xASCII char '65' (= 'A')
        self.assertTrue( np.array_equal(p2dh.parse_calls_into_simplecalls(all_ref_removed_bases.copy(), ref_idx[0], nts_ascii, pileup_line, pileup_lineid), 
                                        np.full(8, nts_ascii[0], dtype=np.int8) )) ## ref 'A' and 'A' called & 4x removed call --> 8xASCII char '65' (= 'A')
        with self.assertRaises(ValueError): ## ref is none --> raise value error
            p2dh.parse_calls_into_simplecalls(all_ref.copy(), None, nts_ascii, pileup_line, pileup_lineid)
    
if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31

