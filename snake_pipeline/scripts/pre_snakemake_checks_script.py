# -*- coding: utf-8 -*-
"""
Runs checks of experiment_info.yaml, collecting errors and eventually outputting a collected list of things to correct/check.

"""

# %%
''' load libraries '''
import os
import sys, argparse
import yaml

cwd = os.getcwd()
sys.path.append(f"{cwd}/scripts")

from gus_helper_functions import *

'''Define functions to be called for error checking'''
def check_samples_csv_header(sample_config_file_path):
    smpl_csv_dict = ['Path','Sample','FileName','Reference','Group','Outgroup']
    correct_header=','.join(smpl_csv_dict)
    with open(sample_config_file_path, 'r') as fid: 
        for line_id, line in enumerate(fid):
            # check if first line is header.
            if line_id == 0:
                if not line.strip('\n')==(correct_header):
                    error_message = f"Header of csv {sample_config_file_path} file is malformed. \n \
                                     Please ensure you use a comma separator and the header is as follows: {correct_header}"
                    return error_message
            else:
                return 0

def check_ref_genome_paths(parsed_config,list_refG):
    path_to_ref_genome_folder = parsed_config['ref_genome_directory']
    errors_list={'Missing file': [],
                 'Missing directory': []}
    # check path exists
    if not os.path.isdir(f'{path_to_ref_genome_folder}'):
        error_message = [f'Path defined in ref_genome_directory field of experiment_info.yaml does not exist! Please check the path {path_to_ref_genome_folder}']
        return error_message
    for ref_genome in set(list_refG):
        if not os.path.isdir(f'{path_to_ref_genome_folder}/{ref_genome}'):
            errors_list['Missing directory'].append(ref_genome)
        elif not os.path.isfile(f'{path_to_ref_genome_folder}/{ref_genome}/genome.fasta') and not os.path.isfile(f'{path_to_ref_genome_folder}/{ref_genome}/genome.fasta.gz'):
            errors_list['Missing file'].append(ref_genome)
    
    # parse errors across possibly multiple reference genomes
    error_message = []
    for ref_genome in errors_list['Missing directory']:
        error_message += [f"Path for reference {ref_genome} does not exist! Please ensure the following path exists {path_to_ref_genome_folder}/{ref_genome}"]
    for ref_genome in errors_list['Missing file']:
        error_message += [f"genome.fasta or genome.fasta.gz file does not exist for reference {ref_genome}! Please ensure this file exists within {path_to_ref_genome_folder}/{ref_genome}"]
    return error_message

def check_ref_genome_main(parsed_config,list_refG):
    snakemake_sections_to_execute = parsed_config['pipeline_specifications']
    if 'mapping' or 'case' in snakemake_sections_to_execute:
        if parsed_config['ref_genome_directory'] == '/PATH/TO/REFERENCE/GENOMES/DIRECTORY':
            return [f"ERROR - pre_snakemake_checks: No reference genome directory specified in experiment_info.yaml. Please provide one, replacing /PATH/TO/REFERENCE/GENOMES/DIRECTORY in the ref_genome_directory field of experiment_info.yaml"]
        return check_ref_genome_paths(parsed_config,list_refG)
            
def check_sample_paths(list_path,list_splID,list_fileN):
    errors_list=[]
    for path,id,filename in zip(list_path,list_splID,list_fileN):
        try:
            findfastqfile(path,id,filename)
        except ValueError as e:
            errors_list.append(e)
    return errors_list

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
                            Runs automated checks at the start of the pipeline execution. \
                            Will output all collected errors to stdout for end-user to fix prior to pipeline execution.

                            #errorfast
                            '''
)
parser.add_argument("-e", dest="experiment_config_yaml", help="experiment_info.yaml file", required=True, action='store')

args = parser.parse_args()

def main():
    # parse files
    with open(args.experiment_config_yaml, "r") as stream:
        try:
            parsed_config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    sample_config_file_path = parsed_config['sample_table']
    # check formatting of samples csv file, fail if necessary:
    if check_samples_csv_header(sample_config_file_path) != 0:
        print(check_samples_csv_header(sample_config_file_path))
        sys.exit(1)

    samples_csv_dict = read_samples_CSV(sample_config_file_path)
    
    # run error checks
    errors = []
    ref_genome_errors = check_ref_genome_main(parsed_config,samples_csv_dict['Reference'])
    if ref_genome_errors:
        errors.append('')
        errors += [f'ERROR - pre_snakemake_checks: Reference genome path parsing failed! Reference genomes referenced in {sample_config_file_path} have the following errors:']
        errors += ['    ' + str(x) for x in ref_genome_errors]

    sample_path_errors = check_sample_paths(samples_csv_dict['Path'],samples_csv_dict['Sample'],samples_csv_dict['FileName'])
    if len(sample_path_errors) > 0:
        errors.append('')
        errors += [f'ERROR - pre_snakemake_checks: File path parsing failed! File paths referenced in {sample_config_file_path} have the following errors:']
        errors += ['    ' + str(x) for x in sample_path_errors]
    if len(errors) > 0:
        for e in errors:
            print(e)
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
   main()
