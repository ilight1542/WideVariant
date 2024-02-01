#!/bin/bash

## This script runs the tests on all generated sample files 
## NOTE: To generate test files run `bash generate_test_cases.sh` in the `testing/` directory 
## To run this script, run `bash run_python_tests.sh` in the `testing/` directory

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
## Run tests
####################

echo -e "\n\n\n######################"
echo "Running boqtie2 QC tests..."
echo -e "######################\n"

python3 python_testing_suite/test_bowtie2qc.py 

echo -e "\n\n\n######################"
echo "Running vcf2quals tests..."
echo -e "######################\n"
python3 python_testing_suite/test_vcf2quals_snakemake.py 
