# Indexes reference genome for bowtie2
rule refGenome_index: 
    input:
        fasta = REF_GENOME_DIRECTORY+"/{reference}/genome.fasta",
    params:
        REF_GENOME_DIRECTORY+"/{reference}/genome_bowtie2",
    output:
        bowtie2idx = REF_GENOME_DIRECTORY+"/{reference}/genome_bowtie2.1.bt2",
    conda:
        "envs/bowtie2.yaml",
    shell:
        "bowtie2-build -q {input.fasta} {params} ;"