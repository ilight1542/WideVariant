# Runs a QC script to summarize results of bowtie2 mapping
rule bowtie2qc:
    input:
        bt2_logs = expand("results/1-mapping/bowtie2/bowtie2_{sampleID}_ref_{reference}.txt",zip,sampleID=samplescsv_dict['Sample'], reference=samplescsv_dict['Reference'] )
        #get_bt2qc_input,
    output:
        alignment_stats = "results/1-mapping/bowtie2_qc/alignment_stats_ref_{reference}.csv",
    params:
        outfile_noextension = "results/1-mapping/bowtie2_qc/alignment_stats_ref_{reference}",
        mappingstats="results/1-mapping/bowtie2",
    conda:
        "../envs/bowtie2qc.yaml",
    shell:
        " python3 {SCRIPTS_DIRECTORY}/bowtie2qc.py -s {spls} -r {wildcards.reference} -m {params.mappingstats} -o {params.outfile_noextension} ;"
