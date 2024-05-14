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

sys.path.append('../scripts/')
import gus_helper_functions as ghf
test_genome_dir_not_findable='./test_data/gus_test_data/non_findable_genome'
test_genome_dir_full='./test_data/gus_test_data/findable_genome'
test_genome_dir_full_gzipped='./test_data/gus_test_data/findable_genome/gzipped'


class TestMyFunction(unittest.TestCase):
    def test_read_fasta(self):
        with self.assertRaises(ValueError):
            ghf.read_fasta(test_genome_dir_not_findable)
        read_in_fasta = ghf.read_fasta(test_genome_dir_full)
        read_in_fasta_gz = ghf.read_fasta(test_genome_dir_full_gzipped)
        self.assertIsInstance(read_in_fasta, SeqIO.FastaIO.FastaIterator)
        self.assertIsInstance(read_in_fasta_gz, SeqIO.FastaIO.FastaIterator)

    def test_nt_order(self):
        self.assertEqual(len(ghf.get_nt_order()),4)

    def test_chrpos2index(self):
        #TODO: add full testing
        #TODO: test error raising
        #TODO: test >2 contigs behaves correctly
        #TODO: test 1 contig behaves correctly
        print('no tests currently implemented')
        pass

    def test_convert_chrpos_to_abspos(self):
        self.assertEqual(ghf.convert_chrpos_to_abspos('scaf1',0,np.array([0,1000]),np.array(['scaf1','scaf2'])),0)
        self.assertEqual(ghf.convert_chrpos_to_abspos('scaf1',500,np.array([0,1000]),np.array(['scaf1','scaf2'])),500)
        self.assertEqual(ghf.convert_chrpos_to_abspos('scaf2',0,np.array([0,1000]),np.array(['scaf1','scaf2'])),1000)
        self.assertEqual(ghf.convert_chrpos_to_abspos('scaf2',500,np.array([0,1000]),np.array(['scaf1','scaf2'])),1500)
        with self.assertRaises(ValueError):
            ghf.convert_chrpos_to_abspos('scaf1',1000,np.array([0,1000]),np.array(['scaf1','scaf2']))
        pass

    def test_round_half_up(self):
        self.assertEqual(ghf.round_half_up(1),1)
        self.assertEqual(ghf.round_half_up(0),0)
        self.assertEqual(ghf.round_half_up(-1),-1)
        self.assertEqual(ghf.round_half_up(-1.5),-1) # next highest integer (from 0.5 cutoff)
        self.assertEqual(ghf.round_half_up(10.5),11)
        self.assertEqual(ghf.round_half_up(11.5),12)
        self.assertEqual(ghf.round_half_up(12.5),13)
        self.assertEqual(ghf.round_half_up(13.5),14)
        self.assertEqual(ghf.round_half_up(14.5),15)
        self.assertEqual(ghf.round_half_up(10.49),10)
        self.assertEqual(ghf.round_half_up(10.49,1),10.5)


    def test_genomestats(self):
        # test ouptut length is as expected
        self.assertEqual(len(ghf.genomestats(test_genome_dir_full)),3)
        self.assertEqual(len(ghf.genomestats(test_genome_dir_full_gzipped)),3)
        # test that starts are correct
        chrstarts_correct=[0,4]
        self.assertEqual(ghf.genomestats(test_genome_dir_full)[0][0],chrstarts_correct[0])
        self.assertEqual(ghf.genomestats(test_genome_dir_full)[0][1],chrstarts_correct[1])
        # test that total length of record is 8
        self.assertEqual(ghf.genomestats(test_genome_dir_full)[1],8)
        # test that scaf names are good
        scafnames_correct=['test_contig_0','test_contig_1']
        self.assertEqual(ghf.genomestats(test_genome_dir_full)[2][0],np.array(scafnames_correct)[0])
        self.assertEqual(ghf.genomestats(test_genome_dir_full)[2][1],np.array(scafnames_correct)[1])

    def test_p2chrpos(self):
        p=np.array([0,1,2,3,4,5,6,7])
        chrstarts_for_p2chrpos_one_chroms=np.array([0])
        chrstarts_for_p2chrpos_two_chroms=np.array([0,4])
        chrstarts_for_p2chrpos_three_chroms=np.array([0,2,4])
        chrstarts_for_p2chrpos_short_contigs=np.array([0,1])
        # test only one contig
        self.assertTrue(np.all(ghf.p2chrpos(p,chrstarts_for_p2chrpos_one_chroms)[:,1]==p))
        # test two ontigs
        self.assertEqual(len(ghf.p2chrpos(p,chrstarts_for_p2chrpos_two_chroms)),8)
        self.assertEqual(np.sum(ghf.p2chrpos(p,chrstarts_for_p2chrpos_two_chroms)[:,0]==0),4)
        self.assertEqual(np.sum(ghf.p2chrpos(p,chrstarts_for_p2chrpos_two_chroms)[:,0]==1),4)
        # test three contigs
        self.assertEqual(len(np.unique(ghf.p2chrpos(p,chrstarts_for_p2chrpos_three_chroms)[:,0])),3)
        # test short contigs
        self.assertEqual(np.sum(ghf.p2chrpos(p,chrstarts_for_p2chrpos_short_contigs)[:,0]==0),1)

if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
