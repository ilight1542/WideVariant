#!/bin/bash
#SBATCH --job-name widevariant.main
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH -o logs/0_widevariant_main_%j.log
#SBATCH -e logs/0_widevariant_main_%j.log
#SBATCH --mail-user=UPDATE_EMAIL_ADDRESS
#SBATCH --mail-type=FAIL

# Activate conda environment (may need to change name of env)
source activate widevariant

snakemake --profile .

# Print "Done!!!" at end of main log file
echo Done!!!
