rule spades:
    input:
        fastq1=rules.sickle2050.output.fq1o,
        fastq2=rules.sickle2050.output.fq2o,
    params:
        outdir="Assembly/spades/{sampleID}"
    conda:
        "../envs/spades_sm.yaml"
    threads: 16
    output:
        fasta="Assembly/spades/{sampleID}/contigs.fasta", # produced by spades''
    shell:
        "spades.py -m 500 -k 21,33,55,77 --phred-offset 33 --careful -t {threads} -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir}"