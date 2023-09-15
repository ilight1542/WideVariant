# Runs a QC script to summarize results of bowtie2 mapping
rule bowtie2qc:
    input:
        get_bt2qc_input,
    output:
        alignment_stats = "1-Mapping/bowtie2_qc/alignment_stats_ref_{reference}.csv",
    params:
        outfile_noextension = "1-Mapping/bowtie2_qc/alignment_stats_ref_{reference}",
    conda:
        "envs/bowtie2qc.yaml",
    shell:
        "python3 {SCRIPTS_DIRECTORY}/bowtie2qc.py -s {spls} -r {wildcards.reference} -d {CURRENT_DIRECTORY} -o {params.outfile_noextension}"