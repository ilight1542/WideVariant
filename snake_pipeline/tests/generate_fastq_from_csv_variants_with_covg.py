# %%
import argparse
import sys
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO,SeqRecord

parser = argparse.ArgumentParser(description='Generator for SNPs from a given CSV. Fields, -i')
parser.add_argument('-i', '--input_variants_csv', metavar='csv', \
    required=True, help="variants csv")
parser.add_argument('-c', '--input_coverage_csv', metavar='csv', \
    required=True, help="covg matrix csv")
parser.add_argument('-r', '--reference', metavar='fasta', \
    required=True, help='fasta to project mutations onto')
parser.add_argument('-l', '--length', metavar='length', \
    required=False, help='length of reads to generate (each end)', \
    default= 150)
parser.add_argument('-I', '--indels_rate', metavar='rate', \
    required=False, help='indel mutation rate', \
    default= 0)
parser.add_argument('-N', '--num_reads', metavar='reads', \
    required=False, help='number of reads to generate per sample', \
    default= 10000)

args=parser.parse_args()

# %%
def add_variants(mutations_to_add,pos_in_genome_of_mutations,refgenome,samplename):
    contig_names, _, _, contig_seqs = parse_ref(refgenome)
    c_name_to_c_seqs=dict(zip(contig_names,contig_seqs))
    contig_pos = mutations_to_add.columns
    for _,sample_variants in mutations_to_add.iterrows():




    for c_name in c_name_to_c_seqs:
        ref_genome_seqs[index_in_record] = str(record.seq)
        ref_genome_chr_names[index_in_record] = [record.id,str(record.description)+' ancestral reconstruction']
    for mutation,pos in zip(mutations_to_add,pos_in_genome_of_mutations):
        ref_genome_seqs[chr]=ref_genome_seqs[chr][:pos] + mutation + ref_genome_seqs[chr][pos+1:]
    seqs_to_output=[]
    for index_in_dicts in range(len(ref_genome_seqs)):
        seq=ref_genome_seqs[index_in_dicts]
        name=ref_genome_chr_names[index_in_dicts][0]
        description=ref_genome_chr_names[index_in_dicts][1]
        seqs_to_output.append(SeqRecord.SeqRecord(Seq(seq), id=name, description=description))
    SeqIO.write(seqs_to_output, f'{samplename}.fasta', "fasta")

def generate_mutated_fastas(raw_variants_csv, reference):
    raw_variants=pd.read_csv(raw_variants_csv,header=1,index=True)
    for samplename,mutations in parsed_csv.iteritems():
        add_variants(mutations,positions,reference,samplename)

def generate_reads_from_fasta(samplename, length, indels_rate, num_reads):
    if indels_rate > 0:
        subprocess.run([f"wgsim -N {num_reads} -r {indels_rate} -R 1 -1 {length} -2 {length} {samplename}.fasta {samplename}_R1.fq {samplename}_R2.fq"],shell=True)
    else: 
        subprocess.run([f"wgsim -N {num_reads} -R 0 -1 {length} -2 {length} {samplename}.fasta {samplename}_R1.fq {samplename}_R2.fq"],shell=True)

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


def generate_reads_per_contig_necessary_for_coverages(input_coverage_csv,lengths,reference):
    parsed_coverages=pd.read_csv(input_coverage_csv,header=1,index=True)
    coverages_for_sample_dict={}
    contig_names,contig_lengths,_,_=parse_ref(reference)
    total_length_readpairs=lengths*2
    for _,r in parsed_coverages.iterrows():
        samplename=r.index
        covg_total=r['Covg total']
        reads_needed_per_contig=[ (r[c_name]*covg_total*c_length)/total_length_readpairs for c_name,c_length in zip(contig_names,contig_lengths) ]
        coverages_for_sample_dict[samplename] = reads_needed_per_contig
    return coverages_for_sample_dict

def main(input_csv,length,reference,indels_rate, num_reads):
    generate_mutated_fastas(parsed_csv, reference)
    for samplename in parsed_csv:
        generate_reads_from_fasta(samplename, length, indels_rate, num_reads)
# %%
if __name__ == '__main__':
    main(input_csv,length,reference,indels_rate)

