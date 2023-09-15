# Turn fastq files into fasta files
rule FQ2FA:
    input:
        fq1o = rules.sickle2050.output.fq1o,
    output:
        fa1o="tmp/{sampleID}_1.fa",
    shell:
        # set +o pipefail; necessary to prevent pipefail (zcat runs but head is done)
        "set +o pipefail; "
        "gzip -cd {input.fq1o} | scripts/fq2fa_sed.sh /dev/stdin > {output.fa1o} ;"


# Run kraken (on forward read file only)
rule kraken2:
    input:
        fa1o = rules.FQ2FA.output.fa1o, # assessment based only on fwd reads
    output:
        kraken_report="Kraken/kraken2/{sampleID}_krakenRep.txt",
        seq_results="Kraken/kraken2/{sampleID}_krakSeq.txt.gz",
    params:
        krakenbracken_db=KRAKEN_BRACKEN_DB, 
    threads: 4,
    conda:
        "envs/kraken_bracken_sm.yaml",
    shell:
        "kraken2 --threads {threads} "
        "--db {params.krakenbracken_db} {input} "
        "--report {output.kraken_report} |gzip > {output.seq_results} "


# Run bracken
rule bracken:
    input:
        kraken_report = rules.kraken2.output.kraken_report,
    output:
        bracken_rep="Kraken/bracken/{sampleID}.bracken",
    params:
        krakenbracken_db=KRAKEN_BRACKEN_DB,
        read_length=KRAKEN_BRACKEN_DB_RL,
    conda:
        "envs/kraken_bracken_sm.yaml",
    shell:
        "scripts/bracken -d {params.krakenbracken_db} -r {params.read_length} -i {input.kraken_report} -o {output.bracken_rep}"