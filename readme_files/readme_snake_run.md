# How to run the Snakemake pipeline


## Outline

This page describes the steps you should take to run the Snakemake pipeline, which processes raw sequencing reads (removes adapters and trims reads based on quality), aligns reads to a reference genome, identifies candidate SNV positions, and saves data arrays with stats on the alignment of each sample at each candidate SNV position. 



## A very brief intro to Snakemake 

If you have not used Snakemake before, please [watch this brief intro video](https://www.youtube.com/watch?v=UOKxta3061g) and refer to the [official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for more information. It will be helpful if you familiarize yourself with the concepts of Snakemake rules, wildcards, and DAGs before proceeding. In brief:
* A *rule* is a small steps that specifies how to create sets of output files from sets of input files. A rule will often employ a specific bioinformatics tool (e.g. bowtie2 for alignment of reads to a reference genome). 
* A *wildcard* is used to generalize a rule to be applicable to a number of datasets. For example, if the same rule will be executed once for each sample, you might use a sample name wildcard when defining input and output files for that rule.
* A *DAG (directed acyclic graph)* refers to the type of flowchart used by Snakemake to determine which rules need to be executed on which files to get from your set of input files to your set of output files. The nodes are rules and the edges are input/output files.


## 0. Install conda (if necessary)

### Installing conda

Conda is a package management system and environment management system. The Snakemake workflow uses conda to manage bioinformatics tools and python packages needed to execute Snakemake rules.

Follow the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for installing Miniconda.


## 1. Download WideVariant and configure

### Copy the template Snakemake directory

All necessary files are available in `snakemake_pipeline/`. The following files you will need to interact with to run the pipeline.
* `snakemake_pipeline/samples.csv` - csv file with one row per sample
* `snakemake_pipeline/experiment_info.yaml` - definded filepaths and variables relevant to this experiment
* `snakemake_pipeline/config.yaml` - information on cluster parameters and what computing resources are required by different workflow steps. Defaults should be reasonable for most slurm systems, but there may be specific modifications that you need to make in order for your HPC and WideVariant to get along.

The following files and directories contain code that should not be modified by the usual end-user. Doing so could break the pipeline or modify its excution in unexpected ways!
* `snakemake_pipeline/Snakefile` - Snakemake workflow
* `snakemake_pipeline/rules/` - rules or small groupings of rules, run within the Snakemake workflow
* `snakemake_pipeline/scripts/` - scripts called by snakemake
* `snakemake_pipeline/envs/` - conda environments used by snakemake

Make a copy of all of these files and directories in the directory on the cluster in which you wish to run the Snakemake workflow. Git clone recommended!

### Modify the necessary input files

#### `snakemake_pipeline/experiment_info.yaml`
NOTE: if planning to run the test data, many of these variables are already defined appropriately.

Mandatory changes:
* change `experiment_name`
* change `sample_table` to path to samples.csv if different than default
* change `ref_genome_directory` to the path which holds your reference genome folders
  * NOTE: Each unique reference genome defined in samples.csv should correspond to a folder name within the `ref_genome_directory` path.
  * For  instance reference genome is `Ypestis_CO92` in samples.csv, and the path defined in `ref_genome_directory` is `/my/path`, the overall structure should be: `/my/path/Ypestis_CO92/` with file `genome.fasta` or `genome.fasta.gz`.
* change `pipeline_specifications` from `["all"]` default, if necessary 
  * This determines which parts of the pipeline snakemake will run
  * `pipeline_specifications` options are: 'mapping', 'case', 'assembly', 'bracken', 'all'
	*  mapping: process reads and align them to a reference genome
    *  case: identify candidate SNVs and generate candidate mutation table
    *  assembly: generate annotated assemblies for each sample
    *  bracken: estimate abundances of taxa in sample
	*  all: run all steps
* if desired, change `generate_normalized_coverage_matrix` and/or `generate_raw_coverage_matrix` to "True"
* change `krakenbracken_db` to the path which holds your kraken bracken DB

#### `snakemake_pipeline/config.yaml`
NOTE: this file may not be renamed

Mandatory changes:
* Put your own email in line `"mail-user":"YOUR_EMAIL_HERE",`, this can also be done at runtime when calling ./WideVariant
* ensure that configfile points to experiment_info.yaml for this experiment (relative path is fine)
* change mkdir, --output and --error to point to log folder for this experiment

Optional changes:
* Uncomment the following lines if you would like to specify cluster partition:
	* `#   --partition={resources.partition}`
	* `#   - partition=""` ; fill in name of partition(s)
* Modify the list of nodes to exclude in line `"exclude":"node327,node088",`
* Change dry-run or unlock parameters if desired 
* Change rule resource requirements


## 2. Modify project-specific files

Some of the files in the template directory must be modified according to the guidelines below.

#### `snakemake_pipeline/samples.csv`

Mandatory changes:
* Fill in the csv with one row per sample and populate all columns
	* Path: the path to the directory where your raw data files are located
		* If you sequenced your sample more than once, then you can list directories to your sequencing data separated by spaces, e.g. "/path/to/seq/run/1/ /path/to/seq/run/2/". Note that the FileName in each directory must be the same.
	* Sample: the name of the sample
		* Constraints: Sample names must be unique
	* FileName: the name of the raw data file for the sample without the file extension and without part of the filename that indicates fwd/rev reads 
		* Example: if your data fwd+rev files are "MySampleBarcode_R2.fq.gz" and "MySampleBarcode_R2.fq.gz", then you should put "MySampleBarcode" as the FileName).
		* Constraints: Raw data files must already be demultiplexed. There must be a file with forward reads and another file with reverse reads (paired end sequencing only).
	* Reference: the name of the reference genome(s) to which Snakemake should align your sample
		* Constraints: The reference genome name must correspond to the name of a directory with a file called "genome.fasta".
	* Group: a string that indicates which samples should be analyzed together
		* The Snakemake workflow will identify candidate SNVs that differentiate samples within the same group. One candidate mutation table will be produced for each group.
		* Constraints: Entry must be a string with no spaces or forbidden characters. All samples in the same group must be aligned to the same reference genome.
	* Outgroup: a boolean value indicating whether or not your sample is an outgroup sample
		* The Snakemake workflow will only identify candidate SNVs that differentiate ingroup samples from each other. It will NOT identify candidate SNVs that differentiate only ingroup samples from outgroup samples.
		* Constraints: Entry must be "0" or "1".



## 3. Run WideVariant on your samples

### Initiate an interactive session on the cluster
This is highly recommended, as the installation of conda environments can be very slow on login nodes. In addition to extra compute power, dispatching the Snakemake from an interactive session will protect against network error and disconnection-related challenges such as corrupted/incomplete conda environments. It may also be recommended to create the initial conda environment prior to calling `./WideVariant`. Commands are below.

Initialize interactive session (Your cluster and/or use case may require different arguments): 
`srun --pty -t 0-10:00:0 -p [PARTITION] -n 1 -c 8 --mem 16G /bin/bash`

Create conda environment:
`conda env create -n widevariant --file envs/widevariant.yaml`

### or Run WideVariant from a compute node or from head node.

Alternatively, the snakemake can be initialized from a non-interactive compute node by setting `run_in_screen` to `False` in experiment_info.yaml. It also can be run directly from a head node if your HPC allows small computations to be carried out on these nodes.

### Run a dry run

A dryrun will catch any issues with incorrect file paths, missing input files, or improper syntax in your Snakefile. A successful dryrun will print out a list of which Snakemake rules will be executed and how many times each one will need to be executed (e.g. the rule `cutadapt` will be executed N times where N is the number of samples you put in `samples.csv`)

You can perform a dryrun by uncommenting `# dryrun: True` in config.yaml before execution. 
Here is an example job tally where there are 8 samples, 4 of which are aligned to one reference genome and the other 4 of which are aligned to another reference genome:

```
Job stats:
job                              count    min threads    max threads
-----------------------------  -------  -------------  -------------
all                                  1              1              1
bowtie2                              8              1              1
bowtie2qc                            2              1              1
candidate_mutation_table             2              1              1
candidate_mutation_table_prep        2              1              1
combine_positions                    2              1              1
combine_positions_prep               2              1              1
cutadapt                             8              1              1
make_data_links                      8              1              1
mpileup2vcf                          8              1              1
pileup2diversity_matrix              8              1              1
refGenome_index                      2              1              1
sam2bam                              8              1              1
sam2bam_cleanup                      8              1              1
samtools_idx                         2              1              1
sickle2050                           8              1              1
variants2positions                   8              1              1
vcf2quals                            8              1              1
total                               95              1              1

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

### Batch the Snakemake workflow to the cluster

Ensure you are in your snakemake_pipeline directory.
1) Start an interactive session if running in screen!
2) `./WideVariant` 

### Track progress

As jobs execute, new directories, `logs/` and `results/` will be made. All final and intermediate results will be in the nested folders in `results/`. 

You can track the workflow's progress in the main log file, `logs/0_widevariant_main.log`. We recommend `watch tail -50 logs/0_widevariant_main.log` if you want to follow along.

### Troubleshooting 

If you have updated your email address either at runtime or within the config.yaml, you will get an email about any failed jobs.

Alternatively, to find jobs that have failed, you can read the main log file `logs/0_widevariant_main.log` and scan for error message, or `grep Error logs/0_widevariant_main.log`. Alternatively, you can print a list of jobs that have failed in the terminal (`sacct | grep FAILED` for failed jobs and `sacct | grep TIMEOUT` for jobs that have timed out). 

To check out the log file for a failed job, type `less logs/*[JOB_ID]*` into the terminal, which will show you the log file contents including any error message. (Type `q` to escape back to the terminal when you are done reading).

If you find that jobs are failing because they ran out of compute resources, then you can change the resource requests in `config.yaml`.

If the main Snakemake workflow is cancelled, you may need to unlock the `Snakefile` to run again. If you do not, there will be an error in `logs/0_widevariant_main.log` indicating that the `Snakefile` is locked. Unlock it with the following:
1) Uncomment `# unlock: True` in `config.yaml`
2) `./WideVariant`
3) Re-comment `# unlock: True` in `config.yaml`

## 4. Download data arrays

Once your Snakemake workflow is completed, you can download the output files for local analysis.

### Output files

Here is a summary of each output file created by this Snakemake workflow:
* `results/2-Case/candidate_mutation_table/group_{groupID}_candidate_mutation_table.npz` 
	* Candidate SNV positions along with alignment stats at each position for each sample
	* One file per Group in `samples.csv`
	* Data objects include:
		* `sample_names` - numpy array of sample names (strings)
			* Dimensions: (number of samples)
		* `p` - numpy array of positions on the genome with candidate SNVs
			* Dimensions: (number of candidate SNV positions)
		* `counts` - numpy array of the number of fwd and rev reads representing each nucleotide at each position for each sample
			* Dimensions: (number of samples) x (number of candidate SNV positions) x 8
		* `quals` - numpy array of the quality scores at each position for each sample
			* Dimensions: (number of samples) x (number of candidate SNV positions)
		* `in_outgroup` - numpy array of boolean representing if each sample was an outgroup sample
			* Dimensions: (number of samples)
		* `indel_counter` - numpy array of statistics on the number of reads supporting insertions and deletions at each position for each sample
			* Dimensions: (number of samples) x (number of candidate SNV positions) x 2
		* `coverage_stats` - (number of samples) x 14 

IMPORTANT NOTE: The above summary represents what these data objects will look like once the output of the Snakemake step is harmonized with the expected inputs of the local analysis step. This has not been updated yet--see [Wishlist for snakemake pipeline upgrades](readme_snake_wishlist.md). The names of the objects are the same, but some of the types and dimensions are different. You can convert them into the right types/dimensions with a quick function at the beginning of the local analysis script. 

If you modified `experiment_info.yaml` to generate coverage matrices, you will get additional output files:
* Please note that raw coverage matrices are often big files.
* `results/2-Case/candidate_mutation_table/group_{cladeID}_coverage_matrix_raw.npz`
	* A coverage matrix with raw read coverage 
	* Dimensions: (number of samples) x (number of positions on genome)
* `results/2-Case/candidate_mutation_table/group_{cladeID}_coverage_matrix_norm.npz`
	* A coverage matrix with normalized read coverage 
	* Dimensions: (number of samples) x (number of positions on genome)



## Alternative Snakemake workflows

The Snakefile provided here can execute alternative workflows, outlined briefly below.

### Assembly workflow

The assembly workflow generates annotated assemblies for each sample.

How to run:
* Change the `pipeline_specifications` variable in `experiment_info.yaml` to include `"assembly"`or `"all"`
* Necessary  fields in `samples.csv`: Path,Sample,FileName,Reference (Reference will not actually be used)

### Taxonomic abundance workflow

The taxonomic abundance workflow generates estimates abundances by taxa in each sample.

How to run:
* Change the `pipeline_specifications` variable in `experiment_info.yaml` to include `"bracken"` or `"all"`
* Necessary  fields in `samples.csv`: Path,Sample,FileName,Reference (Reference will not actually be used)


## Table of contents

[Main Lieberman Lab pipeline README](../README.md)
* [Snakemake pipeline](readme_snake_main.md)
	* [How to run the snakemake pipeline](readme_snake_run.md)
	* [Technical details about the snakemake pipeline](readme_snake_rules.md)
	* [Wishlist for snakemake pipeline upgrades](readme_snake_wishlist.md)
	* [Helpful hints for using the command line](readme_snake_basics.md)
* [Local analysis](readme_local_main.md)
	* [How to run the local analysis script](readme_local_run.md)
	* [Wishlist for local analysis upgrades](readme_local_wishlist.md)
	* [Python best practices](readme_local_best.md)

