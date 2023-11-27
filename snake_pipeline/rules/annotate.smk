# Annotate assembly using bakta
rule bakta:
    input:
        infile=rules.spades.output.fasta,
    params:
        outdir="results/assembly/bakta/{sampleID}",
        db=BAKTA_DATABASE,
    threads: 4
    output:
        stats="results/assembly/bakta/{sampleID}/bakta_out.txt",
        fasta="results/assembly/bakta/{sampleID}/bakta_out.fna",
        faa="results/assembly/bakta/{sampleID}/bakta_out.faa",
        gff="results/assembly/bakta/{sampleID}/bakta_out.gff3",
    conda:
        "../envs/bakta.yaml"
    shell:
        " bakta --db {params.db} --threads {threads} --min-contig-length 500 --output {params.outdir} --prefix bakta_out {input.infile} --force ; "

# Generate annotation stats for each genome
rule bakta_sumstats:
    input:
        fastq=rules.sickle.output.fq1o,
        fasta=rules.bakta.output.fasta,
        assembly_stats=rules.bakta.output.stats,
    output:
        stats="results/assembly/bakta/{sampleID}/sumStats_assembly_annotation.tsv",
    conda:
        "../envs/bakta_sum.yaml",
    shell:
        "python3 {SCRIPTS_DIRECTORY}/bakta_get_sum_stats_script.py -s {wildcards.sampleID} -f {input.fastq} -c {input.fasta} -a {input.assembly_stats} -o {output.stats} "
