#options for testing:
"""
Additional testing step in snakemake pipeline, that runs tests on all outputs
"""

## NOTE: different count matrix sizes have not been tested yet!
## NOTE: some functions are very similar and cound be merged together!
## NOTE: var and unvar pos of ingroup in parse_in_out_variables are generated per testcase but not for individual reference and group --> will fail if multiple refs/groups will be included in test set!

import sys
import unittest
import os
import glob
import Bio.SeqIO as SeqIO
import pandas as pd
import numpy as np


sys.path.append('../scripts/')
import gus_helper_functions as ghf
import build_candidate_mutation_table as bcmtf

maxFQ = -30 
cmt_counts_mat_size = 8
cmt_indelcounter_mat_size = 2
cmt_coverage_stats_shape = 3
cmt_max_cov_to_consider = 200
counts_nt_order = 'ATCGatcg'

class TestMyFunction(unittest.TestCase):
    
    def get_var_pos_raw_dict(self,file, ingroup_samples, experiment_name):
            with open(file, 'r') as fid: 
                header = fid.readline().strip().split(',')[1:] ## read header without sample column
                var_pos_raw_dict_of_exp = {} ## associate for each sample the variant per position
                var_samples_per_loc_dict = {location: 0 for location in header} ## count variants across ingroup samples
                var_pos_ingroup = []
                unvar_pos_ingroup = []
                for line in fid:
                    line = line.strip().split(',')
                    smplname = line[0]
                    smplname_long = f'{experiment_name}_{smplname}'
                    variants_in_sample = line[1:]
                    var_pos_raw_dict_of_exp[smplname_long] = {location: int(variant) for location, variant in zip(header, variants_in_sample)}
                    ## check if sample is ingroup
                    if smplname_long in ingroup_samples:
                        for location, variant in zip(header, variants_in_sample):
                            var_samples_per_loc_dict[location] += int(variant)
                for location in header: 
                    if (var_samples_per_loc_dict[location] == 0) or (var_samples_per_loc_dict[location] == len(ingroup_samples)):
                        unvar_pos_ingroup.append(location) ## get all positions not covered by ingroup to check if those excluded
                    else:
                        var_pos_ingroup.append(location) ## collapse across all ingroup samples if variable positions
            return var_pos_raw_dict_of_exp, var_pos_ingroup, unvar_pos_ingroup
    
    def update_moved_files_in_snakemake_generated_txtfiles(self, infile, outfile, old_path_string, new_path_string):
        with open(infile, 'r') as old_pathnames_file, open(outfile, 'w') as updated_pathnames_out:
            for line in old_pathnames_file:
                # Replace the old path with the new path
                modified_line = line.replace(old_path_string, new_path_string)
                updated_pathnames_out.write(modified_line)
    
    def extract_reference_allele_of_genome(self, var_pos_raw_dict, reference_path):
        first_key = list(var_pos_raw_dict.keys())[0] ## var_pos_raw_dict has samplenames as key and per sample name the same number of locations (keys) and variants (values)
        genome_dict = {}
        variant_ref_dict = {}
        with open(reference_path + 'genome.fasta', 'r') as fasta:
            for record in SeqIO.parse(fasta, 'fasta'):
                genome_dict[record.id] = str(record.seq)
        for location in var_pos_raw_dict[first_key].keys():
            contig, pos = location.split('_')
            if contig in genome_dict.keys():
                variant_ref_dict[location] = genome_dict[contig][int(pos)-1]
            else:
                print('ERROR: contig from variant raw file is not existent in reference genome!')
                sys.exit(0)
        return variant_ref_dict

    def parse_in_out_variables(self):
        cwd=os.getcwd()
        results_subdir = 'test_data/results'
        testcases = glob.glob(f'{cwd}/{results_subdir}/*') #/0-used_input_data/samples.csv')
        arguments_for_testing = {}
        var_pos_raw_dict = {}
        var_pos_ingroup_dict = {}
        unvar_pos_ingroup_dict = {}
        expected_output_dict = {}
        for testcase in testcases:
            experiment_name = testcase.split('/')[-1]
            ## read in variants_raw file used for simulation of SNVs to dictionary
            samplecsv = pd.read_csv(f'{testcase}/0-used_input_data/samples.csv')
            ingroup_samples = samplecsv.loc[samplecsv['Outgroup'] == 0, 'Sample'].values

            ## get positions which are variable in genome
            var_pos_raw_dict[experiment_name], var_pos_ingroup_dict[experiment_name], unvar_pos_ingroup_dict[experiment_name] = self.get_var_pos_raw_dict(f'{testcase}/0-used_input_data/variant_file.csv', ingroup_samples, experiment_name)
            
            ## get unique sets of reference and group
            samplecsv_uniq_ref_grp_pairs = samplecsv[['Reference', 'Group']].drop_duplicates()
            arguments_for_testing[experiment_name] = []
            expected_output_dict[experiment_name] = []
            for reference, groupid in samplecsv_uniq_ref_grp_pairs[['Reference', 'Group']].values:
                samplecsv_of_refgroup = samplecsv[(samplecsv['Reference'] == reference) & (samplecsv['Group'] == groupid)]
                sampleID = samplecsv_of_refgroup["Sample"].values
                outgroup = samplecsv_of_refgroup["Outgroup"].values

                path_to_samplenames_file = f'{testcase}/2-case/temp/group_{groupid}_string_sampleID_names.txt'
                path_to_allpositions = f'{testcase}/2-case/temp/group_{groupid}_allpositions.npz'
                ## modify in file the path names as they have been moved!
                path_to_quals_file = f'{testcase}/2-case/temp/group_{groupid}_string_qual.txt'
                path_to_quals_moved_paths_file = f'{testcase}/2-case/temp/group_{groupid}_string_qual_moved_filenames.txt'
                self.update_moved_files_in_snakemake_generated_txtfiles(path_to_quals_file, path_to_quals_moved_paths_file, 
                                                                        old_path_string = 'results/1-mapping',
                                                                        new_path_string = f'test_data/results/{experiment_name}/1-mapping')
                path_to_diversity_file = f'{testcase}/2-case/temp/group_{groupid}_string_diversity.txt'
                path_to_diversity_moved_paths_file = f'{testcase}/2-case/temp/group_{groupid}_string_diversity_moved_filenames.txt'
                self.update_moved_files_in_snakemake_generated_txtfiles(path_to_diversity_file, path_to_diversity_moved_paths_file, 
                                                                        old_path_string = 'results/1-mapping',
                                                                        new_path_string = f'test_data/results/{experiment_name}/1-mapping')    
                
                path_to_outgroup_file = f'{testcase}/2-case/temp/group_{groupid}_string_outgroup_bool.txt'
                path_to_cmt_out = f'{results_subdir}/{experiment_name}/2-case/candidate_mutation_table/group_{groupid}_candidate_mutation_table.npz' ## need to be relative path!
                path_to_cov_raw_out = f'{results_subdir}/{experiment_name}/2-case/candidate_mutation_table/group_{groupid}_coverage_matrix_raw.npz' ## need to be relative path!
                path_to_cov_norm_out = f'{results_subdir}/{experiment_name}/2-case/candidate_mutation_table/group_{groupid}_coverage_matrix_norm.npz' ## need to be relative path!
                arguments_for_testing[experiment_name].append([path_to_allpositions, path_to_samplenames_file, path_to_outgroup_file, 
                                                                path_to_quals_moved_paths_file, path_to_diversity_moved_paths_file, path_to_cmt_out, 
                                                                path_to_cov_raw_out, path_to_cov_norm_out])
                # arguments_for_testing[experiment_name].append([path_to_allpositions, path_to_samplenames_file, path_to_outgroup_file, 
                #                                                 path_to_quals_moved_paths_file, path_to_diversity_moved_paths_file, path_to_cmt_out, 
                #                                                 path_to_cov_raw_out, path_to_cov_norm_out, flag_cov_raw, flag_cov_norm])

                ## get genome length, num chromosomes and length of each chromosome
                reference_path = f'{cwd}/{reference}/'
                [chr_starts,genome_length,scaf_names] = ghf.genomestats(reference_path)
                ref_nt_on_var_pos = self.extract_reference_allele_of_genome(var_pos_raw_dict[experiment_name], reference_path)
                expected_output_dict[experiment_name].append([sampleID, outgroup, chr_starts, genome_length, scaf_names, ref_nt_on_var_pos])
        return arguments_for_testing, var_pos_raw_dict, var_pos_ingroup_dict, unvar_pos_ingroup_dict, expected_output_dict

    def execute_process_sample_names_tests(self,args_dict_of_testcase,expected_output_of_exp, test_case):
        ## extract sample names from dicts
        path_to_sample_files = args_dict_of_testcase[1]
        expected_sampleids, _, _, _, _, _ = expected_output_of_exp
        
        samplenames_obs, num_samples_obs = bcmtf.process_sample_names(path_to_sample_files)
        ## check if lists have same length and no duplication of sample names occured
        with self.subTest(msg=f'build_CMT_process_sample_names_{test_case}__expected_num_samples'):
            self.assertEqual(num_samples_obs, len(expected_sampleids))
        with self.subTest(msg=f'build_CMT_process_sample_names_{test_case}__no_duplicate_sampleid'):
            self.assertEqual(len(np.unique(samplenames_obs)), num_samples_obs)
        ## check if all expected samples are in samplenames obs
        with self.subTest(msg=f'build_CMT_process_sample_names_{test_case}__all_expected_sampleids_stored'):
            self.assertTrue(all([samplename in samplenames_obs for samplename in expected_sampleids]))
        return num_samples_obs

    def execute_process_positions_tests(self,args_dict_of_testcase,var_pos_ingroup_of_exp, unvar_pos_ingroup_of_exp, expected_output_of_exp, test_case):
        ## extract path to pos file from dicts
        path_to_allpositions = args_dict_of_testcase[0]
        _, _, chr_starts, genome_length, scaf_names, _ = expected_output_of_exp
        positions_obs = bcmtf.process_positions(path_to_allpositions)
        ## check if no position is larger than genome length 
        
        if len(positions_obs) > 0:
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__expected_positions_length'):
                self.assertTrue(max(positions_obs)<=genome_length)
            ## check if no position is smaller than 1 (position are 1-based!)
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__no_pos_smaller_than_1'):
                self.assertTrue(all(positions_obs > 0))
            ## check if any position is duplicated
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__no_duplicated_pos'):
                self.assertEqual(len(np.unique(positions_obs)), len(positions_obs))
            ## convert positions to pos on chromosome 
            chrpos_obs = ghf.p2chrpos(positions_obs, chr_starts)
            chromo_obs = chrpos_obs[:, 0]
            chromo_name_obs = scaf_names[chromo_obs-1]
            pos_obs = chrpos_obs[:, 1]
            chr_pos_obs = [chr_name + '_' + str(pos) for chr_name, pos in zip(chromo_name_obs, pos_obs)]
            ## check if all variable positions of ingroup samples given for data generation are identified
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__var_ingroup_pos'):
                self.assertTrue(all([pos_expected in chr_pos_obs for pos_expected in var_pos_ingroup_of_exp]))
            ## check if all unvariable positions of ingroup samples are not called
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__unvar_ingroup_pos'):
                self.assertTrue(all([pos_expected not in chr_pos_obs for pos_expected in unvar_pos_ingroup_of_exp]))
        else:
            ## check if no variable positions have been expected
            with self.subTest(msg=f'build_CMT_process_positions_{test_case}__no_var_ingroup_pos_expected'):
                self.assertEqual(len(var_pos_ingroup_of_exp), 0)
        return positions_obs

    def execute_process_outgroup_boolean_file_tests(self,args_dict_of_testcase,expected_output_of_exp, test_case):
        ## extract outgroup bools from dicts
        path_to_outgroup_bool_file = args_dict_of_testcase[2]
        _, outgroup, _, _, _, _ = expected_output_of_exp
        
        outgroup_bool_obs = bcmtf.process_outgroup_boolean_file(path_to_outgroup_bool_file)
        ## check if number of samples is correct
        with self.subTest(msg=f'build_CMT_process_outgroup_boolean_file_{test_case}__num_entries'):
            self.assertEqual(outgroup_bool_obs.shape[1], len(outgroup))
        ## check if number of ingroup samples is the expected
        with self.subTest(msg=f'build_CMT_process_outgroup_boolean_file_{test_case}__num_ingroup_samples'):
            self.assertEqual(np.sum(outgroup_bool_obs == False), np.sum(outgroup == 0))
        with self.subTest(msg=f'build_CMT_process_outgroup_boolean_file_{test_case}__num_outgroup_samples'):
            self.assertEqual(np.sum(outgroup_bool_obs), np.sum(outgroup)) ## sum of True and 1 --> no need for conditional check like above

    def execute_process_quals_tests(self,args_dict_of_testcase,num_samples_obs, positions_obs,expected_output_of_exp, test_case):
        ## extract path to quals file from dicts
        path_to_quals_file = args_dict_of_testcase[3]
        sampleID, _, _, _, _, _ = expected_output_of_exp
        quals = bcmtf.process_quals(path_to_quals_file, num_samples_obs, positions_obs)
        ## check if quals matix shape is correct
        with self.subTest(msg=f'build_CMT_process_quals_{test_case}__quals_shape_correct'):
            self.assertEqual(np.shape(quals), (len(positions_obs), len(sampleID))) 
        ## check if all quals are not nan
        with self.subTest(msg=f'build_CMT_process_quals_{test_case}__no_nan_quals'):
            self.assertFalse(np.any(np.isnan(quals))) 
        ## NOTE: quals from paths in path_to_quals_file have been checked to be genomelength long and have quals attribute (see test_vcf2quals)

    def execute_process_counts_tests(self,args_dict_of_testcase,num_samples_obs, positions_obs, expected_output_of_exp, cmt_counts_mat_size, cmt_indelcounter_mat_size, test_case):
        ## extract path to diversity file from dicts
        path_to_diversity_file = args_dict_of_testcase[4]
        sampleID, _, _, genome_length, _, _ = expected_output_of_exp
        counts_obs, indel_counter_obs, all_coverage_per_bp_obs = bcmtf.process_counts(path_to_diversity_file, num_samples_obs, positions_obs)
        ## check if path_to_diversity_file has at least one entry and equals length of samples expected
        with open(path_to_diversity_file, 'r') as f:
            paths_to_diversity_files = f.read().splitlines()
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__path_to_diversity_file_1entry'):
            self.assertTrue( len(paths_to_diversity_files) > 0 )    
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__path_to_diversity_file_numSamplesEntries'):
            self.assertEqual( len(paths_to_diversity_files), len(sampleID) )    

        ## check if counts has correct shape 
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__counts_shape_correct'):
            self.assertEqual( np.shape(counts_obs), (cmt_counts_mat_size, len(positions_obs), len(sampleID)) )
        ## check if no nans are in counts
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__no_nan_counts'):
            self.assertFalse(np.any(np.isnan(counts_obs))) 
        ## check if indel_counter has correct shape 
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__indel_counter_shape_correct'):
            self.assertEqual( np.shape(indel_counter_obs), (cmt_indelcounter_mat_size, len(positions_obs), len(sampleID)) )
        ## check if no nans are in indel_counter
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__no_nan_indel_counter'):
            self.assertFalse(np.any(np.isnan(indel_counter_obs))) 
        ## check if all_coverage_per_bp has correct shape 
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__all_coverage_per_bp_shape_correct'):
            self.assertEqual( np.shape(all_coverage_per_bp_obs), (genome_length, len(sampleID)) )
        ## check if no nans are in all_coverage_per_bp
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__no_nan_all_coverage_per_bp'):
            self.assertFalse(np.any(np.isnan(all_coverage_per_bp_obs))) 
        ## NOTE: Output depends on wgsim output and correct position association was checked in pileup2 diversity test
        return all_coverage_per_bp_obs
    
    def execute_calculate_coverage_distribution_exceed_maxcov_tests(self):
        
        max_coverage_for_test = 10
        all_coverage_per_bp = np.array([[ 8,  0], [15,  3], [11,  3], [ 9,  5], [ 8,  4], [ 7,  1], [10, 14], [ 0,  7], [14, 12], [ 7, 11]]) ## 2 samples with 10 nt long genome
        num_samples = 2
        coverage_stats_exp = np.array([[8.5, 8.9, 3.961060464067672], [4.5, 6.0, 4.58257569495584]], dtype='uint')
        coverage_distribution_exp = np.array([[1,0,0,0,0,0,0,2,2,1,4], [1,1,0,2,1,1,0,1,0,0,3]], dtype='uint')
        coverage_stats_obs, coverage_distribution_obs = bcmtf.calculate_coverage_distribution(all_coverage_per_bp, num_samples, max_coverage_for_test)
        ## check if coveragestats is returning correct coverage_distribution
        with self.subTest(msg=f'build_CMT_process_counts__coverage_stats_above_max_cov'):
            self.assertTrue( np.array_equal(coverage_stats_obs, coverage_stats_exp) )
        ## check if coveragestats is returning correct coverage_distribution
        with self.subTest(msg=f'build_CMT_process_counts__coverage_dist_above_max'):
            self.assertTrue( np.array_equal(coverage_distribution_obs, coverage_distribution_exp) )
        

    def execute_calculate_coverage_distribution_tests(self,all_coverage_per_bp_obs, num_samples_obs, expected_output_of_exp, cmt_coverage_stats_shape, cmt_max_cov_to_consider, test_case):

        sampleID, _, _, _, _, _ = expected_output_of_exp
        coverage_stats_obs, coverage_distribution_obs = bcmtf.calculate_coverage_distribution(all_coverage_per_bp_obs, num_samples_obs)
        ## select the dimension of the coverage distribution matrix (either the highest coverage if below max_threshold, or max_threshold!)
        max_cov_across_samples = np.max( all_coverage_per_bp_obs )
        max_dimension_cov_dist = np.min( [max_cov_across_samples, cmt_max_cov_to_consider])
        ## check if coverage_stats_obs has correct shape 
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__coverage_stats_shape_correct'):
            self.assertEqual( np.shape(coverage_stats_obs), (len(sampleID), cmt_coverage_stats_shape) )
        ## check if no nans are in coverage_stats_obs
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__no_nan_all_coverage_stats'):
            self.assertFalse(np.any(np.isnan(coverage_stats_obs))) 
        ## check if coverage_distribution_obs has correct len
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__coverage_distribution_len_correct'):
            self.assertEqual( np.shape(coverage_distribution_obs)[0], len(sampleID) )
        ## check if coverage_distribution_obs is not longer than max given coverage
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__coverage_distribution_max_cov_correct'):
            self.assertEqual( np.shape(coverage_distribution_obs)[1],  max_dimension_cov_dist+1 )
        ## check if no nans are in coverage_distribution_obs
        with self.subTest(msg=f'build_CMT_process_counts_{test_case}__no_nan_all_coverage_per_bp'):
            self.assertFalse(np.any(np.isnan(coverage_distribution_obs))) 

    def execute_normalize_coverage_matrix_invar_samples_tests(self):
        all_coverage_per_bp_1smpl_invar = np.array([[ 3,  0], [ 3,  3], [ 3,  3], [ 3,  5], [ 3,  4]]) ## 2 samples with genomelength 5; 1 sample invariable coverage
        all_coverage_per_bp_all_invar = np.array([[ 3,  3], [ 3,  3], [ 3,  3], [ 3,  3], [ 3,  3]]) ## 2 samples with genomelength 5; both invariable
        all_coverage_per_bp_pos_invar = np.array([[ 5,  5], [ 3,  3], [ 8,  8], [10, 10], [ 7,  7]]) ## 2 samples with genomelength 5; invariable per position
        array_cov_norm_scaled_1smpl_invar_exp = np.array([[ 1000, -1000], [ 0, 0], [ 0, 0], [-1000,  1000], [-1000,  1000]])
        array_cov_norm_scaled_all_invar_exp = np.zeros((5, 2), dtype = 'uint') ## no variance within samples --> first normalization sets to 0 --> all values are 0!
        array_cov_norm_scaled_pos_invar_exp = np.zeros((5, 2), dtype = 'uint') ## no variance within samples --> first normalization sets to 0 --> all values are 0!

        array_cov_norm_scaled_1smpl_invar_obs = bcmtf.normalize_coverage_matrix(all_coverage_per_bp_1smpl_invar)
        array_cov_norm_scaled_all_invar_obs = bcmtf.normalize_coverage_matrix(all_coverage_per_bp_all_invar)
        array_cov_norm_scaled_pos_invar_obs = bcmtf.normalize_coverage_matrix(all_coverage_per_bp_pos_invar)
        ## check if array_cov_norm_scaled_1smpl_invar_obs is as expected
        with self.subTest(msg=f'build_CMT_normalize_coverage_matrix__expected_output_1sample_invariable'):
            self.assertTrue( np.array_equal(array_cov_norm_scaled_1smpl_invar_obs, array_cov_norm_scaled_1smpl_invar_exp) )
        ## check if array_cov_norm_scaled_all_invar_exp is all 0 
        with self.subTest(msg=f'build_CMT_normalize_coverage_matrix__expected_output_1sample_invariable'):
            self.assertTrue( np.array_equal(array_cov_norm_scaled_all_invar_obs, array_cov_norm_scaled_all_invar_exp) )
        ## check if array_cov_norm_scaled_pos_invar_obs is all 0
        with self.subTest(msg=f'build_CMT_normalize_coverage_matrix__expected_output_1sample_invariable'):
            self.assertTrue( np.array_equal(array_cov_norm_scaled_pos_invar_obs, array_cov_norm_scaled_pos_invar_exp) )


    def execute_normalize_coverage_matrix_tests(self,all_coverage_per_bp_obs, test_case):
        
        array_cov_norm_scaled_obs = bcmtf.normalize_coverage_matrix(all_coverage_per_bp_obs)
        ## check if array_cov_norm_scaled and all_coverage_per_bp_obs have same shape
        with self.subTest(msg=f'build_CMT_nomalize_coverage_matrix_{test_case}__expected_shape_double_norm_cov_mat'):
            self.assertEqual( np.shape(array_cov_norm_scaled_obs), np.shape(all_coverage_per_bp_obs) )
        ## check if array_cov_norm_scaled contains no nans
        with self.subTest(msg=f'build_CMT_nomalize_coverage_matrix_{test_case}__no_nan_double_norm_cov_mat'):
            self.assertFalse(np.any(np.isnan(array_cov_norm_scaled_obs))) 
        ## NOTE: Due to rounding errors one cannot check if matrix sums up to 0 (or if the np.std is 0!)
        

    def check_correct_npzfile_existence(self, path, npzattribute, shape, testmessage_prefix, return_mat = False):
        ## check if file exists
        with self.subTest(msg=f'{testmessage_prefix}_file_exists'):
            self.assertTrue(os.path.exists(path))
        npz_file = np.load(path)
        npz_file_attribute = npz_file.files
        ## check if attribute is present as attribute in file
        with self.subTest(msg=f'{testmessage_prefix}_file_expected_attribute'):
            self.assertTrue(npzattribute in npz_file_attribute)
        npz_file_mat = npz_file[npzattribute]
        ## check if npzfile has expected shape
        with self.subTest(msg=f'{testmessage_prefix}_expected_shape'):
            self.assertTrue(np.shape(npz_file_mat), shape)
        if return_mat: 
            return npz_file_mat

    def check_CMT_file_existance(self, path_to_candidate_mutation_table, expected_attribute_in_final_cmt, testmessage_prefix, return_cmt_file = True):
        with self.subTest(msg=f'{testmessage_prefix}__CMT_file_exists'):
            self.assertTrue(os.path.exists(path_to_candidate_mutation_table))
        cmt_file = np.load(path_to_candidate_mutation_table)
        cmt_file_attributes = cmt_file.files
        for attribute in expected_attribute_in_final_cmt:
            ## check if expected attributes are present as attribute in cmt
            with self.subTest(msg=f'{testmessage_prefix}__CMT_file_expected_attribute_{attribute}'):
                self.assertTrue(attribute in cmt_file_attributes)
        if return_cmt_file:
            return cmt_file

    def check_CMT_sample_names(self, sample_names_obs, shape, samplenames_exp, testmessage_prefix):
        ## check if cmt_matrices have expected shapes
        with self.subTest(msg=f'{testmessage_prefix}_expected_num_samples'):
            self.assertEqual(np.shape(sample_names_obs), shape)
        ## check if all expected samples are in samplenames obs
        with self.subTest(msg=f'{testmessage_prefix}_all_expected_sampleids_stored'):
            self.assertTrue(all([samplename in sample_names_obs for samplename in samplenames_exp]))

    def check_CMT_positions(self, positions_obs, chr_pos_obs, var_pos_ingroup_exp, unvar_pos_ingroup_exp, testmessage_prefix):
        ## check if cmt_matrices have at least number of positions as ingroup has variable positions
        with self.subTest(msg=f'{testmessage_prefix}_min_num_of_positions'):
            self.assertTrue(np.shape(positions_obs)[0] >= len(var_pos_ingroup_exp))
        ## check if all variable positions are given
        with self.subTest(msg=f'{testmessage_prefix}_all_var_ingroup_pos_encoded'):
            self.assertTrue( all([var_position_exp in chr_pos_obs for var_position_exp in var_pos_ingroup_exp]) )
        ## check if none of unvariable positions are given
        with self.subTest(msg=f'{testmessage_prefix}_no_unvar_ingroup_pos_encoded'):
            self.assertTrue( all([unvar_position_exp not in chr_pos_obs for unvar_position_exp in unvar_pos_ingroup_exp]) )

    def get_maNT_from_counts(self, counts_obs, counts_nt_order):
        ## extract the column idx of counts matrix to sum up per nt
        counts_nt_order_upper_list = np.array([nt.upper() for nt in counts_nt_order]) 
        nt_idx_dict = {nt_upper: np.where(nt_upper == counts_nt_order_upper_list) for nt_upper in np.unique(counts_nt_order_upper_list)}

        ## get the major allele called in dataset
        maNT_array = np.empty(np.shape(counts_obs)[:2], dtype='object') ## initialize array with same shape as counts but no nt resolution (2D instead of 3D!)
        for i in range(counts_obs.shape[0]):
            for j in range(counts_obs.shape[1]):
                max_sum = 0
                max_group = ''
                # Calculate the sum for each nt and find the max
                for group, indices in nt_idx_dict.items():
                    group_sum = counts_obs[i, j, indices[0]].sum()
                    if group_sum > max_sum:
                        max_sum = group_sum
                        max_group = group
                # Store the group with the maximum sum
                maNT_array[i, j] = max_group
        return maNT_array

    def check_CMT_counts(self, counts_obs, sample_names_obs, chr_pos_obs, shape, var_pos_raw_dict, variant_ref_dict, counts_nt_order, testmessage_prefix):
        ## check if cmt_matrix has correct shape
        with self.subTest(msg=f'{testmessage_prefix}__expected_shape'):
            self.assertEqual(np.shape(counts_obs), shape)
        ## check if no nans are in cmt_matrix counts
        with self.subTest(msg=f'{testmessage_prefix}__no_nan'):
            self.assertFalse(np.any(np.isnan(counts_obs)))
        
        if shape[1] > 0: ## at least one observed position
            maNT_array = self.get_maNT_from_counts(counts_obs, counts_nt_order)
            ## check if at variant position, ref is not major allele and an major allele is given!
            for samplename, variants_on_sample in var_pos_raw_dict.items():
                ## get sample idx of counts matrix
                sample_idx = np.where(samplename == sample_names_obs)[0][0]
                for location, is_variant in variants_on_sample.items():
                    ## get position idx of counts matrix
                    pos_idx = np.where(location == chr_pos_obs)[0][0]
                    ## get reference as expected
                    ref_on_pos = variant_ref_dict[location]
                    ## get major allele nt:
                    if is_variant == 1:
                        ## variant --> check if reference nt does not align with called allele 
                        with self.subTest(msg=f'{testmessage_prefix}__var_pos_maNT_unequal_refNT'):
                            self.assertTrue(maNT_array[sample_idx, pos_idx] != ref_on_pos)
                    else:
                        ## no variant --> check if reference nt _does_ align with called allele 
                        with self.subTest(msg=f'{testmessage_prefix}__unvar_pos_maNT_is_refNT'):
                            self.assertEqual(maNT_array[sample_idx, pos_idx], ref_on_pos)
                

    def check_CMT_quals(self, quals_obs, sample_names_obs, chr_pos_obs, shape, var_pos_raw_dict, maxFQ, testmessage_prefix):
        ## check if cmt_matrix has correct shape
        with self.subTest(msg=f'{testmessage_prefix}_expected_shape'):
            self.assertEqual(np.shape(quals_obs), shape)
        ## check if no nans are in cmt_matrix counts
        with self.subTest(msg=f'{testmessage_prefix}_no_nan'):
            self.assertFalse(np.any(np.isnan(quals_obs)))
        ## check if at given position for given sample with variant the qual score is below threshold 
        for sampleid, variants_dict in var_pos_raw_dict.items():
            sample_idx = np.where(sampleid == sample_names_obs)[0]
            quals_of_sample_obs = quals_obs[sample_idx, :].flatten()
            for quals_of_sample_pos_obs, pos_obs in zip(quals_of_sample_obs, chr_pos_obs):
                if pos_obs in variants_dict.keys():
                    if variants_dict[pos_obs] == 1: ## check if variant 
                        ## check if variant has FQ score lower than threshold FQ score
                        with self.subTest(msg=f'{testmessage_prefix}_variant_lower_than_thresh_FQscore'):
                            self.assertTrue(quals_of_sample_pos_obs < maxFQ)

    def check_CMT_outgroup(self, outgroup_obs, sample_names_obs, shape, samplenames_exp, outgroups_exp, testmessage_prefix):
        ## check if cmt_matrix has expected shapes
        with self.subTest(msg=f'{testmessage_prefix}_expected_shape'):
            self.assertEqual(np.shape(sample_names_obs), shape)
        ## check if all sampleids have correct outgroup association
        for samplename_exp, outgroup_exp in zip(samplenames_exp, outgroups_exp):
            sample_idx = np.where(sample_names_obs == samplename_exp)[0]
            outgroup_of_smpl_obs = int(outgroup_obs[sample_idx][0])
            with self.subTest(msg=f'{testmessage_prefix}_sampleid_correct_outgroup_association'):
                self.assertEqual(outgroup_of_smpl_obs, outgroup_exp)
    
    def check_CMT_check_shape_nans(self, mat, shape, testmessage_prefix):
        ## check if cmt_matrix has expected shape
        with self.subTest(msg=f'{testmessage_prefix}_expected_shape'):
            self.assertEqual(np.shape(mat), shape)
        ## check if no nan value existent
        with self.subTest(msg=f'{testmessage_prefix}__no_nan'):
            self.assertFalse(np.any(np.isnan(mat)))
    
    def check_CMT_coverage_dist_obs(self, mat, shape, testmessage_prefix):
        ## check if cmt_matrix has expected shape
        with self.subTest(msg=f'{testmessage_prefix}_expected_sample_dimension_length'):
            self.assertEqual(np.shape(mat)[0], shape[0])
        with self.subTest(msg=f'{testmessage_prefix}_expected_cov_dimension_length'):
            self.assertTrue(np.shape(mat)[1] <= shape[1])
        ## check if no nan value existent
        with self.subTest(msg=f'{testmessage_prefix}__no_nan'):
            self.assertFalse(np.any(np.isnan(mat)))

    def execute_build_CMT_tests(self,args_dict_of_testcase, cmt_counts_mat_size, cmt_indelcounter_mat_size, cmt_max_cov_to_consider, var_pos_raw_dict, var_pos_ingroup_of_exp, unvar_pos_ingroup_of_exp, maxFQ, expected_output_of_exp, test_case):
        ## NOTE: Some redundant testing for shape of files and the correct encoding of positions done
        ## but necessary to ensure that everything as expected is encoded in CMT in the end!!
        expected_attribute_in_final_cmt = ['sample_names', 'p', 'counts', 'quals', 'in_outgroup', 'indel_counter', 'coverage_stats', 'coverage_dist']
        sampleID, outgroup, chr_starts, genome_length, scaf_names, variant_ref_dict = expected_output_of_exp
        [_,_,_,_,_, path_to_candidate_mutation_table, path_to_cov_mat_raw, path_to_cov_mat_norm] = args_dict_of_testcase
        
        ## check output for the different possible flags
        for cov_raw_flag in [True, False]:
            if cov_raw_flag:
                message_string_raw = '_covraw'
            else:
                message_string_raw = ''
            for cov_norm_flag in [True, False]:
                bcmtf.build_CMT(*args_dict_of_testcase, cov_raw_flag, cov_norm_flag) 
                if cov_norm_flag:
                    message_string_comb = message_string_raw + '_covnorm'
                else:
                    message_string_comb = message_string_raw
                ####
                ## raw coverage matrix
                if cov_raw_flag:
                    cov_raw = self.check_correct_npzfile_existence(path_to_cov_mat_raw, 'cov_mat_raw', (len(sampleID), genome_length), f'build_CMT_build_CMT_{test_case}{message_string_comb}__covraw', return_mat=True)
                    max_cov_obs = np.max(cov_raw)
                    max_shape_cov_dist = int(np.min( [max_cov_obs, cmt_max_cov_to_consider] )+1)
                else:
                    max_shape_cov_dist = int(cmt_max_cov_to_consider+1)
                ####
                ## norm coverage matrix
                if cov_norm_flag:
                    cov_norm = self.check_correct_npzfile_existence(path_to_cov_mat_norm, 'cov_mat_norm_scaled', (len(sampleID), genome_length), f'build_CMT_build_CMT_{test_case}{message_string_comb}__covnorm', return_mat=True)
                    
                ####
                ## general CMT 
                ## check CMT for correct dimensions and existence of file
                cmt_file = self.check_CMT_file_existance(path_to_candidate_mutation_table, expected_attribute_in_final_cmt, f'build_CMT_build_CMT_{test_case}{message_string_comb}', return_cmt_file = True)
                
                ####
                ## Sample names in CMT
                sample_names_obs = cmt_file['sample_names']
                self.check_CMT_sample_names(sample_names_obs, (len(sampleID), ), sampleID, f'build_CMT_build_CMT_{test_case}{message_string_comb}__sample_names_CMT')
                
                ####
                ## positions in CMT
                positions_obs = cmt_file['p']
                ## convert positions to pos on chromosome 
                chrpos_obs = ghf.p2chrpos(positions_obs, chr_starts)
                chromo_name_obs = scaf_names[chrpos_obs[:, 0]-1]
                pos_obs = chrpos_obs[:, 1]
                chr_pos_obs = np.array([chr_name + '_' + str(pos) for chr_name, pos in zip(chromo_name_obs, pos_obs)])
                self.check_CMT_positions(positions_obs, chr_pos_obs, var_pos_ingroup_of_exp, unvar_pos_ingroup_of_exp, f'build_CMT_build_CMT_{test_case}{message_string_comb}__positions_CMT')
                num_obs_positions = len(positions_obs) 
                

                ## counts in CMT
                counts_obs = cmt_file['counts']
                self.check_CMT_counts(counts_obs, sample_names_obs, chr_pos_obs, (len(sampleID), num_obs_positions, cmt_counts_mat_size), var_pos_raw_dict, variant_ref_dict, counts_nt_order, f'build_CMT_build_CMT_{test_case}{message_string_comb}__counts_CMT')
    
                ## quals in CMT
                quals_obs = cmt_file['quals']
                self.check_CMT_quals(quals_obs, sample_names_obs, chr_pos_obs, (len(sampleID), num_obs_positions), var_pos_raw_dict, maxFQ, f'build_CMT_build_CMT_{test_case}{message_string_comb}__quals_CMT') ## Note sample_names_obs, chr_pos_obs needed to identify correct idx of quals matrix 
                
                ## outgroup in CMT
                outgroup_obs = cmt_file['in_outgroup']
                self.check_CMT_outgroup(outgroup_obs, sample_names_obs, (len(sampleID), ), sampleID, outgroup, f'build_CMT_build_CMT_{test_case}{message_string_comb}__outgroup_CMT')
    
                ## indel counter in CMT
                indel_counter_obs = cmt_file['indel_counter']
                self.check_CMT_check_shape_nans(indel_counter_obs, (len(sampleID), num_obs_positions, cmt_indelcounter_mat_size), f'build_CMT_build_CMT_{test_case}{message_string_comb}__indel_counter_CMT')

                ## coverage stats in CMT
                coverage_stats_obs = cmt_file['coverage_stats']
                self.check_CMT_check_shape_nans(coverage_stats_obs, (len(sampleID), cmt_coverage_stats_shape), f'build_CMT_build_CMT_{test_case}{message_string_comb}__coverage_stats_CMT')

                ## coverage distribution in CMT
                coverage_dist_obs = cmt_file['coverage_dist']
                if cov_raw_flag:
                    ## max_shape_cov_dist reflects actual bound of cov dimension 
                    self.check_CMT_check_shape_nans(coverage_dist_obs, (len(sampleID), max_shape_cov_dist), f'build_CMT_build_CMT_{test_case}{message_string_comb}__coverage_dist_CMT')
                else:
                    ## max_shape_cov_dist reflects just upper bound --> just see if coverage dimension is in bounds of <= max_shape_cov_dist
                    self.check_CMT_coverage_dist_obs(coverage_dist_obs, (len(sampleID), max_shape_cov_dist), f'build_CMT_build_CMT_{test_case}{message_string_comb}__coverage_dist_CMT')

    def test_outputs(self):
        # run a separate instance of test_outputs for each test dataset (and each refgenome within those test datasets)
        args_dict, var_pos_raw_dict, var_pos_ingroup_dict, unvar_pos_ingroup_dict, expected_output_dict = self.parse_in_out_variables()

        self.execute_calculate_coverage_distribution_exceed_maxcov_tests()
        self.execute_normalize_coverage_matrix_invar_samples_tests()
        for test_case in args_dict.keys():
            var_pos_raw_of_testcase = var_pos_raw_dict[test_case] 
            var_pos_ingroup_of_testcase = var_pos_ingroup_dict[test_case] 
            unvar_pos_ingroup_of_testcase = unvar_pos_ingroup_dict[test_case]
            ## loop over the different group/references in the arguments
            for single_args_dict, single_expected_output_dict in zip(args_dict[test_case], expected_output_dict[test_case]):
                ## test each function of combine positions script independently 
                num_samples_obs = self.execute_process_sample_names_tests(single_args_dict, single_expected_output_dict, test_case) ## num samples observed needed for later evaluation (used as input of subsequent functions)!
                positions_obs = self.execute_process_positions_tests(single_args_dict, var_pos_ingroup_of_testcase, unvar_pos_ingroup_of_testcase, single_expected_output_dict, test_case) ## positions observed needed for later evaluation (used as input of subsequent functions)!
                self.execute_process_outgroup_boolean_file_tests(single_args_dict, single_expected_output_dict, test_case)
                self.execute_process_quals_tests(single_args_dict, num_samples_obs, positions_obs, single_expected_output_dict, test_case) 
                all_coverage_per_bp_obs = self.execute_process_counts_tests(single_args_dict, num_samples_obs, positions_obs, single_expected_output_dict, cmt_counts_mat_size, cmt_indelcounter_mat_size, test_case) ## all coverage bp observed needed for later evaluation (used as input of subsequent functions)!
                self.execute_calculate_coverage_distribution_tests(all_coverage_per_bp_obs, num_samples_obs, single_expected_output_dict, cmt_coverage_stats_shape, cmt_max_cov_to_consider, test_case)
                self.execute_normalize_coverage_matrix_tests(all_coverage_per_bp_obs, test_case)
                self.execute_build_CMT_tests(single_args_dict, cmt_counts_mat_size, cmt_indelcounter_mat_size, cmt_max_cov_to_consider, 
                                            var_pos_raw_of_testcase, var_pos_ingroup_of_testcase, unvar_pos_ingroup_of_testcase, maxFQ, single_expected_output_dict, test_case)
            print(f'build_CMT oython script checked on {test_case}\n\n\n')   



if __name__ == '__main__':
    unittest.main() ## when run as script
else:
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31unittest.main(argv=['first-arg-is-ignored'], exit=False) ## when run in ipyhton or ijupyter it is needed to ignore return value from sys.argv (first argument), see https://medium.com/@vladbezden/using-python-unittest-in-ipython-or-jupyter-732448724e31
