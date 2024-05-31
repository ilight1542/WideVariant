#!/bin/bash

## This script runs the tests on all generated sample files 
## NOTE: To generate test files run `bash generate_test_cases.sh` in the `testing/` directory 
## To run this script, run `bash run_python_tests.sh` in the `testing/` directory

####################
## FUNCTIONS
####################

## function to invoke test, capture the std out and err and print last three lines (if test was ok)
run_test() {
    local test_script_basename=$1
    local path_to_log=$2
    local timestamp=$3

    ## run test
    {
        python3 python_testing_suite/${test_script_basename}.py 
    } >> "${path_to_log}/${timestamp}_${test_script_basename}.log" 2>&1

    ####
    ## Error reporting
    ## echo last three lines of logfile (time and if test was ok!)
    unittest_report=$(grep -A2 -E "Ran .* test in .*s" "${path_to_log}/${timestamp}_${test_script_basename}.log")
    echo "$unittest_report"

    ## Check if error occured which was not spotted by unittest
    if ( grep -i -E 'error|fail|traceback' ${path_to_log}/${timestamp}_${test_script_basename}.log | \
            grep -v 'File exists' | \
            grep -vq 'Falling back to greedy solver' ) ; then     
    
        echo "An Error occured during testing. Please check the log file (${path_to_log}/${timestamp}_${test_script_basename}.log) for more information "; 
    fi
}

####################
## Activate Conda
####################

## get conda env and activate
conda_envs_path=$(echo $CONDA_EXE | sed 's#bin/conda#envs#g' )
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')

if [[ -d "${conda_envs_path}/widevariant_snakemake_testing" ]];then
    source ${activate_conda_path} widevariant_snakemake_testing
else
    conda env create -n widevariant_snakemake_testing --file testing/widevariant_snakemake_testing.yaml
    source ${activate_conda_path} widevariant_snakemake_testing
fi

####################
## Set global variables
####################

log_file_path='test_data/logs/testing'
timestamp=$(date +"%Y_%m_%d_%H_%M")

####################
## Run tests
####################

echo
echo ${timestamp}
echo "Start testing..."
mkdir -p ${log_file_path}

echo -e "\n\n\n######################"
echo "Running general_snakemake_helpers tests..."
echo -e "######################\n"
script_basename='test_general_snakemake_helpers'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running bowtie2 QC tests..."
echo -e "######################\n"
script_basename='test_bowtie2qc'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running combinepositions tests..."
echo -e "######################\n"
script_basename='test_combinepositions'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running variants2positions tests..."
echo -e "######################\n"
script_basename='test_variants2positions'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running vcf2quals tests..."
echo -e "######################\n"
script_basename='test_vcf2quals_snakemake'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running pileup2diversity_helper_functions tests..."
echo -e "######################\n"
script_basename='test_pileup2diversity_helper_functions'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo -e "\n\n\n######################"
echo "Running build_candidate_mutation_table tests..."
echo -e "######################\n"
script_basename='test_build_candidate_mutation_table'
run_test ${script_basename} ${log_file_path} ${timestamp}
echo -e "${script_basename} done\n"

echo "Tests done"
echo "If some errors have been reported, please check the log files ${log_file_path} and look for 'ERROR', 'FAILED' and 'Traceback'"
echo