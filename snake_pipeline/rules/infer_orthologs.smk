# Get two-column (caldeID,path2faa) input file for ortholog_inference script
rule build_annotation_orthologs_input:
    input:
        prokka_faa=expand("results/assembly/prokka/{sampleID}/prokka_out.faa",sampleID=samplescsv_dict['Sample']),
    params:
        clade_identifier=expand("{sampleID}",sampleID=samplescsv_dict['Sample']),
    output:
        "results/assembly/orthologinfo_filtered/input_files.tsv",
    group:
        "infer_orthologs"
    shell:
        """
        paste <(echo {params.clade_identifier} | {SCRIPTS_DIRECTORY}/sed_nl.sh ) <(echo {input.prokka_faa} | {SCRIPTS_DIRECTORY}/sed_nl.sh ) > {output}
        """

# Infer ortholog info for each identified gene (based on AA sequence) across all clades using CD-HIT
rule run_annotation_orthologs_inference:
    input:
        rules.build_annotation_orthologs_input.output
    params:
        percent_identity="0.9", # percent identity for clustering
        cdhit_mem="8000", # max mem available for cdhit
        output_folder="results/assembly/orthologinfo_filtered/"
    output:
        "results/assembly/orthologinfo_filtered/annotation_orthologs.tsv"
    conda:
        "../envs/cdhit.yaml"
    group:
        "infer_orthologs"
    shell:
        "python3 {SCRIPTS_DIRECTORY}/annotation_orthologs_inference.py -f {input} -p {params.percent_identity} -m {params.cdhit_mem} -o {params.output_folder}"
