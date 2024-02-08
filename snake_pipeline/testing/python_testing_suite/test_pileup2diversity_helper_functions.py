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

sys.path.append('./scripts/')
import gus_helper_functions as ghf


class TestMyFunction(unittest.TestCase):
    def test_generate_nts_ascii(self):
        nt_order=ghf.get_nt_order()
        ord_conversion={A: 65, T: 84, C: 67, G: 71, a: 97, t: 116, c: 99, g: 103}
        self.assertEqual(p2dh.generate_nts_ascii()[:4], [ord_conversion(x) for x in ['A','T','C','G']])
        self.assertEqual(p2dh.generate_nts_ascii()[4:], [ord_conversion(x) for x in ['a','t','c','g']])
    
    def test_clean_calls_from_start_and_end(self):
        no_start_or_end=np.fromstring('....,,,,', dtype=np.int8)
        all_start=np.fromstring('^X.^X.^X.^X.^X,^X,^X,^X,', dtype=np.int8)
        all_end=np.fromstring('$X.$X.$X.$X.$X,$X,$X,$X,', dtype=np.int8)
        mixed_start_end=np.fromstring('$X.^X.$X.^X.$X,^X,$X,^X,', dtype=np.int8)
        # test that lengths are correct for clean calls
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_start),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_end),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(mixed_start_end),8)

    def test_parse_indels_into_data(self):
        #TODO
        test_data=np.array()
        no_start_or_end=np.fromstring('....,,,,', dtype=np.int8)
        all_insertions=np.fromstring('^X.^X.^X.^X.^X,^X,^X,^X,', dtype=np.int8)
        all_deletions=np.fromstring('$X.$X.$X.$X.$X,$X,$X,$X,', dtype=np.int8)
        mixed_indels=np.fromstring('$X.^X.$X.^X.$X,^X,$X,^X,', dtype=np.int8)
        # test that lengths are correct for clean calls
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_start),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_end),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(mixed_start_end),8)

    def parse_calls_into_simplecalls(self):
        #TODO
        no_start_or_end=np.fromstring('....,,,,', dtype=np.int8)
        all_insertions=np.fromstring('^X.^X.^X.^X.^X,^X,^X,^X,', dtype=np.int8)
        all_deletions=np.fromstring('$X.$X.$X.$X.$X,$X,$X,$X,', dtype=np.int8)
        mixed_indels=np.fromstring('$X.^X.$X.^X.$X,^X,$X,^X,', dtype=np.int8)
        # test that lengths are correct for clean calls
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_start),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(all_end),8)
        self.assertEqual(p2dh.clean_calls_from_start_and_end(mixed_start_end),8)
    
if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31

