import argparse
import sys
import pandas as pd
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Generator for SNPs from a given CSV. Fields, -i')
parser.add_argument('-i', '--input', metavar='csv', \
    required=True, nargs=1, help="Lenght of probes to create")
parser.add_argument('-r', '--reference', metavar='fasta', \
    required=True, nargs=1, help='fasta to project mutations onto')
parser.add_argument('-l', '--length', metavar='length', \
    required=False, nargs=1, help='length of reads to generate (each end)', \
    default= 150,)
parser.add_argument('-I', '--indels_rate', metavar='rate', \
    required=False, nargs=1, help='indel mutation rate', \
    default= 0,)
parser.add_argument('-N', '--num_reads', metavar='reads', \
    required=False, nargs=1, help='number of reads to generate per sample', \
    default= 10000,)

args=parser.parse_args()

def add_variants(mutations_to_add,pos_in_genome_of_mutations,refgenome,samplename)
    refgenome = SeqIO.parse(refgenome,'fasta')
    ref_genome_seqs={}
    ref_genome_chr_names={}
    for index_in_record,record in enumerate(refgenome):
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

def generate_mutated_fastas(parsed_csv, reference):
    positions=parsed_csv.index
    for samplename,mutations in parsed_csv.iteritems():
        add_variants(mutations,positions,reference,samplename)


def generate_reads_from_fasta(samplename, length, indels_rate, num_reads):
    if indels_rate > 0:
        subprocess.run([f"wgsim -N {num_reads} -r {indels_rate} -R 1 -1 {length} -2 {length} {samplename}.fasta {samplename}_R1.fq {samplename}_R2.fq"],shell=True)
    else: 
        subprocess.run([f"wgsim -N {num_reads} -R 0 -1 {length} -2 {length} {samplename}.fasta {samplename}_R1.fq {samplename}_R2.fq"],shell=True)

def main(input_csv,length,reference,indels_rate, num_reads):
    parsed_csv=pd.read_csv(input_csv,header=1,index=True)
    generate_mutated_fastas(parsed_csv, reference)
    for samplename in parsed_csv:
        generate_reads_from_fasta(samplename, length, indels_rate, num_reads)


if __name__ == '__main__':
    main(input_csv,length,reference,indels_rate)

