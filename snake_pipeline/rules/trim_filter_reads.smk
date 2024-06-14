# Prepare filtered, clean FASTQ samples
rule cutadapt:
    ## NOTE: for proper identification this must be turned on, but currently causes issues with service output // pipe links between rules and inappropriate resource usages
    # input:
    #     samplecsv = "results/data/{sampleID}/sample_info.csv",
    params:
        fq1 = "results/data/{sampleID}/R1.fq.gz",
        fq2 = "results/data/{sampleID}/R2.fq.gz",
        adapter_illumina_nextera="CTGTCTCTTAT",
        adapter_illumina_universal="AGATCGGAAGAGC",
    output:
        fq1o = "results/preprocessing/{sampleID}/R1_trim.fq.gz",
        fq2o = "results/preprocessing/{sampleID}/R2_trim.fq.gz",
        log = "logs/preprocessing_{sampleID}.txt", # necessary for bowtie2qc
    threads: 4,
    conda:
        "../envs/cutadapt.yaml",
    shell:
        # NEEDS TO BE OF FORMAT FILENAMER1, FILENAMER2, NOITHING ELSE
        "cutadapt -a {params.adapter_illumina_nextera} "
                "-a {params.adapter_illumina_universal} "
                "-A {params.adapter_illumina_nextera} "
                "-A {params.adapter_illumina_universal} "
                "-m 1 --cores={threads} "
                "-o {output.fq1o} -p {output.fq2o} "
                "{params.fq1} {params.fq2} "
                "1> {output.log};"

rule sickle:
    input:
        fq1i = rules.cutadapt.output.fq1o,
        fq2i = rules.cutadapt.output.fq2o,
    params:
        qual=20, # Threshold for trimming based on average quality in a window
        readlen=50, # Threshold to keep a read based on length after trimming
    conda:
        "../envs/sickle-trim.yaml",
    output:
        fq1o = "results/preprocessing/{sampleID}/R1_filt.fq.gz",
        fq2o = "results/preprocessing/{sampleID}/R2_filt.fq.gz",
        fqSo = "results/preprocessing/{sampleID}/filt_sgls.fq.gz",
        log = "logs/preprocessing_{sampleID}.txt", # necessary for bowtie2qc
    shell:
        "sickle pe -f {input.fq1i} -r {input.fq2i} "
                    "-t sanger "
                    "-o {output.fq1o} -p {output.fq2o} "
                    "-s {output.fqSo} "
                    "-g -q {params.qual} -l {params.readlen} -x -n "
                    "1>> {output.log} ;"
        "echo 'sickle command line params: minqual={params.qual} minreadlen={params.readlen}' 1>> {output.log} ;"
