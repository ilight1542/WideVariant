#!/bin/bash

####################
## FUNCTIONS
####################

## function to generate data
generate_data() {
    local experiment_name=$1
    local variant_file=$2
    local coverage_file=$3
    local outgroup_ids=$4
    local genome=${5:-'testing/Smallgenome/genome.fasta'}
    local read_length=${6:-75}

    rm -rf test_data/input_data/${experiment_name}

    python3 ${generate_test_data_pythonscript} \
        -n ${experiment_name} \
        -i ${variant_file} \
        -c ${coverage_file} \
        -o ${outgroup_ids} \
        -r ${genome} \
        -l ${read_length}
}
## function to run snakemake
run_snakemake() {
    ## get into parent folder to envoke snakemake run
    cd ../
    
    rm -rf results/
    rm -rf logs/*
    
    snakemake \
        --cores 1 \
        --use-conda \
        --conda-prefix ${conda_envs_path} \
        --conda-frontend mamba \
        --latency-wait 30 \
        --keep-going \
        --configfile ./experiment_info.yaml \
        --rerun-triggers mtime 
    ## change back to testing folder
    cd testing
}
## function to move data to proper subfolder
move_and_link_data() {
    local experiment_name=$1
    local variant_file=$2
    local coverage_file=$3
    local outgroup_ids=$4
    local genome=${5:-'testing/Smallgenome/genome.fasta'}

    ## move results to correct subfolder
    rm -rf test_data/results/${experiment_name}
    mkdir -p test_data/results/${experiment_name}/0-used_input_data
    mv ../results/* test_data/results/${experiment_name}
    ## generate softlinks from input data for testing 
    ln -s $(pwd)/${variant_file} test_data/results/${experiment_name}/0-used_input_data/variant_file.csv
    ln -s $(pwd)/${coverage_file} test_data/results/${experiment_name}/0-used_input_data/coverage_file.csv
    ln -s $(pwd)/${outgroup_ids} test_data/results/${experiment_name}/0-used_input_data/outgroup_ids.csv
    ln -s $(pwd)/${genome} test_data/results/${experiment_name}/0-used_input_data/genome.fasta
    ln -s $(pwd)/test_data/sample_csvs/${experiment_name}.csv test_data/results/${experiment_name}/0-used_input_data/samples.csv
}

check_for_errors_in_log() {
    local logfile=$1
    
    ## search in case insensitive manner (-i), with regex (-E), exclude certain patterns (-v) and just report if true (-q)
    ## "File exists" warning comes from ln if soft link is already present
    ## "Falling back to greedy solver" warning from snakemake for scheduling
    if ( grep -i -E 'error|fail|traceback' ${logfile} | \
            grep -v 'File exists' | \
            grep -vq 'Falling back to greedy solver' ) ; then     

        echo "An Error occured during data generation. Please check the log file (${logfile}) for more information "; 
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

## run snakemake interactively on test data
generate_test_data_pythonscript='python_testing_suite/generate_test_data.py'
genome='Smallgenome/genome.fasta'
read_length=75
log_file_path='test_data/logs/data_generation'
timestamp=$(date +"%Y_%m_%d_%H_%M")

####################
## Test cases
####################

echo
echo ${timestamp}
echo "Start generating test data..."

mkdir -p ${log_file_path}

####
# Default: ground_truth (no outgroup sample) 
experiment_name='ground_truth'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"

####
# Single sample input (CURRENTLY BREAKS SOME SNAKEMAKE RULES) 
experiment_name='single_sample'
variant_file='data_generation/variants_raw_single_sample.csv'
coverage_file='data_generation/coverage_single_sample.csv'
outgroup_ids='data_generation/outgroup_ids_single_sample.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"

####
# Single outgroup sample (first in csv) 
experiment_name='single_outgroup_first'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_first_outgroup.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"


####
# Single outgroup sample (last in csv) 
experiment_name='single_outgroup_last'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_last_outgroup.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"


####
# Single ingroup sample 
experiment_name='single_ingroup'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_one_ingroup.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"


####
# 2 input samples, one outgroup sample one ingroup sample 
experiment_name='two_input_out_and_in'
variant_file='data_generation/variants_raw_two_samples.csv'
coverage_file='data_generation/coverage_two_samples.csv'
outgroup_ids='data_generation/outgroup_ids_two_samples_one_outgroup_one_ingroup.csv'

{
    generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length}
    run_snakemake 
    move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
} >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
echo -e "${experiment_name} done\n"


####
# Multiple outgroups samples



echo "Data generation done"
echo "If some errors have been reported, please check the log files ${log_file_path} and look for 'ERROR', 'FAILED' and 'Traceback'"
echo