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
    local additional_flag_for_generate_test_data_script=$7

    rm -rf test_data/input_data/${experiment_name}

    python3 ${generate_test_data_pythonscript} \
        -n ${experiment_name} \
        -i ${variant_file} \
        -c ${coverage_file} \
        -o ${outgroup_ids} \
        -r ${genome} \
        -l ${read_length} ${additional_flag_for_generate_test_data_script}
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

run_data_generation_with_logging() {
    local experiment_name=$1
    local variant_file=$2
    local coverage_file=$3
    local outgroup_ids=$4
    local genome=$5
    local read_length=$6
    local log_file_path=$7
    local timestamp=$8
    local additional_flag_for_generate_test_data_script=$9
    
    {
        generate_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${additional_flag_for_generate_test_data_script}
        run_snakemake 
        move_and_link_data ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome}
    } >> "${log_file_path}/${timestamp}_${experiment_name}.log" 2>&1

    check_for_errors_in_log "${log_file_path}/${timestamp}_${experiment_name}.log"
    echo -e "${experiment_name} done\n"
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

# List of test cases
# 1. Default: ground_truth (no outgroup sample) 
# 2. Single sample input (CURRENTLY BREAKS SOME SNAKEMAKE RULES)
# 3. Single outgroup sample (first in csv) 
# 4. Single outgroup sample (last in csv) 
# 5. Single ingroup sample & multiple outgroup samples
# 6. 2 input samples, one outgroup sample one ingroup sample 
# 7. Variants only on outgroup samples
# 8. Not any variant on any sample
# 9. One sample with no variant, one sample only with variants
# 10. One position with only variant, one position with variant in ingroup samples only
# 11. One sample with only one variant which is not shared with any other sample
# 12. Variants on every first and last position of all contigs


echo
echo ${timestamp}
echo "Start generating test data..."

mkdir -p ${log_file_path}

# 12. Variants on every first and last position of all contigs
experiment_name='contig_edges_vars'
variant_file='data_generation/variants_raw_contig_edges.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'
set_contigs_to_be_covered='-e'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp} ${set_contigs_to_be_covered}



####
# 1. Default: ground_truth (no outgroup sample) 
experiment_name='ground_truth'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 2. Single sample input (CURRENTLY BREAKS SOME SNAKEMAKE RULES) 
experiment_name='single_sample'
variant_file='data_generation/variants_raw_single_sample.csv'
coverage_file='data_generation/coverage_single_sample.csv'
outgroup_ids='data_generation/outgroup_ids_single_sample.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 3. Single outgroup sample (first in csv) 
experiment_name='single_outgroup_first'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_first_outgroup.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 4. Single outgroup sample (last in csv) 
experiment_name='single_outgroup_last'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_last_outgroup.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 5. Single ingroup sample & multiple outgroup samples
experiment_name='single_ingroup'
variant_file='data_generation/variants_raw.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_one_ingroup.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 6. 2 input samples, one outgroup sample one ingroup sample 
experiment_name='two_input_out_and_in'
variant_file='data_generation/variants_raw_two_samples.csv'
coverage_file='data_generation/coverage_two_samples.csv'
outgroup_ids='data_generation/outgroup_ids_two_samples_one_outgroup_one_ingroup.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 7. Variants only on outgroup samples
experiment_name='var_only_in_outgroup'
variant_file='data_generation/variants_raw_only_outgroup_var.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_two_outgroups.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 8. Not any variant on any sample
experiment_name='no_var_any_sample'
variant_file='data_generation/variants_raw_not_any_var.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 9. One sample with no variant, one sample only with variants
experiment_name='two_samples_all_var_no_var'
variant_file='data_generation/variants_raw_no_var_all_var.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 10. One position with only variant, one position with variant in ingroup samples only
experiment_name='two_pos_all_var_ingroup_var'
variant_file='data_generation/variants_raw_invariant_position.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids_last_outgroup.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 11. One sample with only one variant which is not shared with any other sample (sample B position contig1 1000)
experiment_name='one_var_uniq_across_samples'
variant_file='data_generation/variants_raw_one_sample_unique_var.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}

####
# 12. Variants on every first and last position of all contigs
experiment_name='contig_edges_vars'
variant_file='data_generation/variants_raw_contig_edges.csv'
coverage_file='data_generation/coverage.csv'
outgroup_ids='data_generation/outgroup_ids.csv'

run_data_generation_with_logging ${experiment_name} ${variant_file} ${coverage_file} ${outgroup_ids} ${genome} ${read_length} ${log_file_path} ${timestamp}



echo "Data generation done"
echo "If some errors have been reported, please check the log files ${log_file_path} and look for 'ERROR', 'FAILED' and 'Traceback'"
echo