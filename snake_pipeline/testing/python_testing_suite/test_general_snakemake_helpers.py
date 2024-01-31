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

test_genome_dir_not_findable='testing/test_data/gus_test_data/non_findable_genome'
test_genome_dir_full='testing/test_data/gus_test_data/findable_genome'
test_genome_dir_full_gzipped='testing/test_data/gus_test_data/findable_genome/gzipped'


class TestMyFunction(unittest.TestCase):
    def test_read_fasta(self):
        self.assertRaises(gus.read_fasta(test_genome_dir_not_findable))
        read_in_fasta = gus.read_fasta(test_genome_dir_full)
        read_in_fasta_gz = gus.read_fasta(test_genome_dir_full_gzipped)
        self.assertIsInstance(read_in_fasta, SeqIO.Seq)
        self.assertIsInstance(read_in_fasta_gz, SeqIO.Seq)

    def test_genomestats(self):
        # test ouptut length is as expected
        self.assertEqual(len(gus.genomestats(test_genome_dir_full)),4)
        self.assertEqual(len(gus.genomestats(test_genome_dir_full_gzipped)),4)
        # test that starts are correct
        chrstarts_correct=[0,4]
        self.assertEqual(len(gus.genomestats(test_genome_dir_full)[1]),np.array())
        # test that total length of record is 8
        self.assertEqual(len(gus.genomestats(test_genome_dir_full)[1]),8)
        # test that scaf names are good
        scafnames_correct=['test_contig_0','test_contig_1']
        self.assertEqual(gus.genomestats(test_genome_dir_full)[2]),np.array(scafnames_correct)

    def test_p2chrpos(self):
        p=np.array([0,1,2,3,4,5,6,7])
        chrstarts_for_p2chrpos_one_chroms=np.array([0])
        chrstarts_for_p2chrpos_two_chroms=np.array([0,4])
        # test only one contig
        self.assertEqual(gus.p2chrpos(p,chrstarts_for_p2chrpos_one_chroms)[:,1],p)
        # test contigs
        self.assertEqual(len(gus.p2chrpos(p,chrstarts_for_p2chrpos_two_chroms)[1,:]),4)

if __name__,'__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
