import argparse
import os, shutil
import sys
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO,SeqRecord,Seq

parser = argparse.ArgumentParser(description='Generator for SNPs from a given CSV. Fields, -i')
parser.add_argument('-n','--name', metavar='str', required=True, help='name of run/experiment')
parser.add_argument('-i', '--input_variants_csv', metavar='csv', \
    required=True, help="variants boolean csv")
parser.add_argument('-c', '--input_coverage_csv', metavar='csv', \
    required=True, help="covg matrix csv")
parser.add_argument('-r', '--reference', metavar='fasta', \
    required=True, help='fasta to project mutations onto')
parser.add_argument('-l', '--length', metavar='length', \
    required=False, help='length of reads to generate (each end)', \
    default= 150)
parser.add_argument('-b', '--basecalls_csv',metavar='csv',required=False,help='specific ATCG mutations for the variant bool positions csv', default=None)
parser.add_argument('-o', '--outgroup_csv',metavar='csv',required=False,help='specific ATCG mutations for the variant bool positions csv', default=None)
parser.add_argument('-e', '--cover_contig_edges', required=False, help='Adjust wgsim generated FastQ files to cover contig edges with mean contig coverage', action='store_true')
parser.add_argument('-s', '--save_tmp', required=False, help='Save temp file directory within tests (tmp)', action='store_true')

args=parser.parse_args()

# %%
def iter_variant(ref):
    nts=np.array(['A','T','C','G'])
    ref_idx=np.where(nts == ref)[0][0]
    return np.take(nts,ref_idx+1,mode='wrap')

def generate_fasta_name(experiment_name, samplename, contig_name):
    return f'{experiment_name}_sample_{samplename}_contig_{contig_name}'

def parse_ref(reference):
    ref_genome = SeqIO.parse( reference, 'fasta' )
    genome_length = 0 # init
    contig_lengths = [] # init
    contig_names = [] # init
    contig_seqs = [] # init

    for record in ref_genome: # loop through contigs
        contig_names.append(record.id)
        genome_length = genome_length + len(record)
        contig_lengths.append(len(record))
        contig_seqs.append(str(record.seq))

    # Turn into numpy arrays
    genome_length = np.asarray( genome_length, dtype=np.int_ )
    contig_names = np.asarray( contig_names, dtype=object )
    contig_lengths = np.asarray( contig_lengths, dtype=object )
    contig_seqs = np.asarray( contig_seqs, dtype=object )

    return [ contig_names, contig_lengths, genome_length, contig_seqs ]

def generate_mutated_fastas(experiment_name, variants_boolean_csv,refgenome,basecalls_csv=None):
    """
    Parses variant positions in the variants boolean csv and outputs a new fasta file for each sample x contig for fasta generation
    """
    if not os.path.isdir("tmp"): 
        os.makedirs("tmp") 
    parsed_variants=pd.read_csv(variants_boolean_csv,header=0,index_col=0)
    if basecalls_csv:
        basecalls = pd.read_csv(variants_boolean_csv,header=0,index_col=0).values
    else: basecalls = None
    contig_pos=[x.split('_') for x in parsed_variants.columns]
    for sample_index,sample_variants in parsed_variants.iterrows():
        contig_names, _, _, contig_seqs = parse_ref(refgenome)
        c_name_to_c_seqs=dict(zip(contig_names,contig_seqs))
        sample_variants_bool = np.array(sample_variants)
        for variant_position in np.where(sample_variants_bool)[0]:
            c_name,c_seq_idx=contig_pos[variant_position]
            c_seq_idx = int(c_seq_idx) - 1 ## -1 for 0-based correction
            ref = c_name_to_c_seqs[c_name][c_seq_idx]
            if basecalls == None:
                alt = iter_variant(ref)
            else:
                alt = basecalls[sample_index,variant_position]
            c_name_to_c_seqs[c_name] = c_name_to_c_seqs[c_name][:c_seq_idx] + alt + c_name_to_c_seqs[c_name][c_seq_idx+1:]
        for c_name in c_name_to_c_seqs:
            seq=c_name_to_c_seqs[c_name]
            seq_to_output = SeqRecord.SeqRecord(Seq.Seq(seq), id=f'{c_name}', description=f'Mutated contig_{c_name}')
            fasta_name=generate_fasta_name(experiment_name, sample_index, c_name)
            SeqIO.write(seq_to_output, f'tmp/{fasta_name}.fasta', "fasta")

def generate_reads_per_contig_necessary_for_coverages(experiment_name, input_coverage_csv,lengths,reference):
    parsed_coverages=pd.read_csv(input_coverage_csv,header=0,index_col=0)
    num_reads_for_sample_dict={}
    coverage_for_sample_dict={}
    contig_names,contig_lengths,_,_=parse_ref(reference)
    total_length_readpairs=int(lengths)*2
    for samplename,r in parsed_coverages.iterrows():
        covg_total=r['Covg total']
        coverage_for_sample = [ (float(r[c_name])*float(covg_total)) for c_name,c_length in zip(contig_names,contig_lengths) ]
        reads_needed_per_contig=[ (float(r[c_name])*float(covg_total)*int(c_length))/total_length_readpairs for c_name,c_length in zip(contig_names,contig_lengths) ]
        for index,num_reads in enumerate(reads_needed_per_contig):
            fasta_name=generate_fasta_name(experiment_name, samplename, contig_names[index])
            num_reads_for_sample_dict[f'{fasta_name}'] = int(num_reads)
            coverage_for_sample_dict[f'{fasta_name}'] = coverage_for_sample[index]
    return num_reads_for_sample_dict, coverage_for_sample_dict

def generate_reads_from_fasta(fasta_file, fastq_r1, fastq_r2, length, num_reads):
    subprocess.run([f"wgsim -N {num_reads} -r 0 -e 0 -R 0 -1 {length} -2 {length} {fasta_file} {fastq_r1} {fastq_r2}"],shell=True)

def combine_reads_across_contigs(experiment_name,sample_names,contig_names):
    if not os.path.isdir(f"test_data/input_data/{experiment_name}"): 
        os.makedirs(f"test_data/input_data/{experiment_name}") 
    for sample in sample_names:
        input_files_R1=[f'tmp/{experiment_name}_sample_{sample}_contig_{c}_R1.fq' for c in contig_names]
        input_files_R2=[f'tmp/{experiment_name}_sample_{sample}_contig_{c}_R2.fq' for c in contig_names]
        command_cat_R1 = f"cat {' '.join(input_files_R1)} > test_data/input_data/{experiment_name}/{experiment_name}_{sample}_R1.fastq"
        command_cat_R2 = f"cat {' '.join(input_files_R2)} > test_data/input_data/{experiment_name}/{experiment_name}_{sample}_R2.fastq"
        # Execute the shell command
        subprocess.run(command_cat_R1, shell=True)
        subprocess.run(command_cat_R2, shell=True)
        subprocess.run(f'gzip test_data/input_data/{experiment_name}/{experiment_name}_{sample}_R1.fastq', shell=True)
        subprocess.run(f'gzip test_data/input_data/{experiment_name}/{experiment_name}_{sample}_R2.fastq', shell=True)

def prepare_test_run_widevariant_call(experiment_name,sample_names,reference,outgroups):
    cwd=os.getcwd()
    reference_full_path=f'{cwd}/{reference}'
    reference_name=reference.split('/')[-2]
    reference_upstream_path=reference_full_path.split(f'/{reference_name}')[0]
    # create new samples.csv for this experiment in tests
    if not os.path.isdir(f"test_data/sample_csvs"):
        os.makedirs(f"test_data/sample_csvs")
    with open(f'test_data/sample_csvs/{experiment_name}.csv','w') as f:
        f.write('Path,Sample,FileName,Reference,Group,Outgroup\n')
        for sample,outgroup in zip(sample_names,outgroups):
            f.write(f'{cwd}/test_data/input_data/{experiment_name},{experiment_name}_{sample},{experiment_name}_{sample},{reference_name},1,{outgroup}\n')
    # modify experiment_info.yaml 
    lines_to_output=[]
    with open('../experiment_info.yaml') as f:
        for l in f:
            l=l.strip()
            if l.startswith('experiment_name'):
                l = f"experiment_name: {experiment_name}"
            if l.startswith('sample_table'):
                l = f'sample_table: testing/test_data/sample_csvs/{experiment_name}.csv'
            if l.startswith('ref_genome_directory'):
                l = f'ref_genome_directory: {reference_upstream_path}'
            lines_to_output.append(l)
    with open('../experiment_info.yaml', 'w') as f:
        for l in lines_to_output:
            f.write(l)
            f.write('\n')

def get_sample_names_contig_names_outgroups(coverage_csv,outgroup_csv):
    parsed_covg=pd.read_csv(coverage_csv,header=0,index_col=0)
    parsed_covg=parsed_covg[parsed_covg.columns[parsed_covg.columns != 'Covg total']]
    sample_names=parsed_covg.index
    contig_names=parsed_covg.columns
    if outgroup_csv == None:
        outgroup_ids = [0 for s in sample_names]
    else:
        parsed_outgroup=pd.read_csv(outgroup_csv,header=0,index_col=0)
        outgroup_ids=parsed_outgroup['Outgroup']
    return sample_names,contig_names,outgroup_ids

def get_contig_edge_sequence(fasta_file):
    with open(fasta_file) as fasta: 
        for record in SeqIO.parse(fasta, "fasta"):
            return record.seq ## single fasta file

def modify_fastq(fastq_file, random_read_id_contig_start, random_read_id_contig_end, contig_seq, rev_complement = False):
    contig_length = len(contig_seq)
    mod_records = []
    with open(fastq_file, 'r') as fastq: 
        for read_id, record in enumerate(SeqIO.parse(fastq, "fastq")):
            if read_id in random_read_id_contig_start:
                start, end = record.id.split('_')[1:3]
                fragment_size = abs(int(end) - int(start))
                len_read = len(record.seq)
                record.id = record.id.replace(f'_{start}_{end}_', f'_1_{fragment_size+1}_')
                record.description = ''
                if not rev_complement:
                    record.seq = contig_seq[:len_read]
                else:
                    record.seq = Seq.reverse_complement(contig_seq[fragment_size-len_read:fragment_size])
            elif read_id in random_read_id_contig_end:
                start, end = record.id.split('_')[1:3]
                fragment_size = abs(int(end) - int(start))
                len_read = len(record.seq)
                record.id = record.id.replace(f'_{start}_{end}_', f'_{contig_length-fragment_size}_{contig_length}_')
                record.description = ''
                if not rev_complement:
                    record.seq = contig_seq[-fragment_size:-fragment_size+len_read]
                else:
                    record.seq = Seq.reverse_complement(contig_seq[-len_read:])
            mod_records.append(record)
    ## overwrite old fastq file with adjusted reads
    with open(fastq_file, 'w') as fastqout: 
        SeqIO.write(mod_records, fastqout, "fastq")

def generate_reads_on_contig_edges(R1_fastq_file, R2_fastq_file, fasta_file, reads_on_contig, num_reads_on_contig_edge):
    if num_reads_on_contig_edge < 2: ## enforce at least 2, if function is called
        num_reads_on_contig_edge = 2
    random_read_id_contig_start, random_read_id_contig_end = np.array_split(np.random.choice(reads_on_contig, int(num_reads_on_contig_edge), replace = False), 2)
    ## get fasta sequence 
    contig_seq = get_contig_edge_sequence(fasta_file)
    ## modify the randomly selected reads to sit on contig edges
    modify_fastq(R1_fastq_file, random_read_id_contig_start, random_read_id_contig_end, contig_seq, rev_complement = True)
    modify_fastq(R2_fastq_file, random_read_id_contig_start, random_read_id_contig_end, contig_seq)
    #print(f'{int(num_reads_on_contig_edge)} reads have been modified to sit at contig edges')

def clean_tmp():
    shutil.rmtree('tmp')

def main(input_experiment_name,input_variants_csv,input_coverage_csv,input_basecalls_csv,input_outgroup_csv,input_length,input_reference,cover_contig_edges,input_save_tmp):
    sample_names,contig_names,outgroup_ids = get_sample_names_contig_names_outgroups(input_coverage_csv,input_outgroup_csv)
    generate_mutated_fastas(input_experiment_name, input_variants_csv, input_reference,input_basecalls_csv)
    num_reads_for_sample_dict, coverage_for_sample_dict = generate_reads_per_contig_necessary_for_coverages(input_experiment_name, input_coverage_csv,input_length,input_reference)
    for sample_contig_to_output in num_reads_for_sample_dict:
        num_reads=num_reads_for_sample_dict[sample_contig_to_output]
        fasta_file = f'tmp/{sample_contig_to_output}.fasta'
        fastq_r1 = f'tmp/{sample_contig_to_output}_R1.fq'
        fastq_r2 = f'tmp/{sample_contig_to_output}_R2.fq'
        generate_reads_from_fasta(fasta_file, fastq_r1, fastq_r2, input_length, num_reads)
        if cover_contig_edges:
            coverage=coverage_for_sample_dict[sample_contig_to_output]
            generate_reads_on_contig_edges(fastq_r1, fastq_r2, fasta_file, num_reads, coverage)
    combine_reads_across_contigs(input_experiment_name,sample_names,contig_names)
    prepare_test_run_widevariant_call(input_experiment_name,sample_names,input_reference,outgroup_ids)
    if not input_save_tmp:
        clean_tmp()

if __name__ == '__main__':
    main(args.name,args.input_variants_csv,args.input_coverage_csv,args.basecalls_csv,args.outgroup_csv,args.length,args.reference,args.cover_contig_edges,args.save_tmp)
