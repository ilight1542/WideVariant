# Get two-column (caldeID,path2faa) input file for ortholog_inference script
rule build_annotation_orthologs_input:
    input:
        prokka_faa=expand("Assembly/prokka/{sampleID}/prokka_out.faa",sampleID=SAMPLE_ls_long),
    params:
        clade_identifier=expand("{sampleID}",sampleID=SAMPLE_ls_long),
    output:
        "Assembly/orthologinfo_filtered/input_files.tsv",
    shell:
        """
        paste <(echo {params.clade_identifier} | scripts/sed_nl.sh ) <(echo {input.prokka_faa} | scripts/sed_nl.sh ) > {output}
        """

# Infer ortholog info for each identified gene (based on AA sequence) across all clades using CD-HIT
rule infer_orthologs:
    input:
        rules.build_annotation_orthologs_input.output
    params:
        percent_identity="0.9", # percent identity for clustering
        cdhit_mem="8000", # max mem available for cdhit
        output_folder="Assembly/orthologinfo_filtered/"
    output:
        "Assembly/orthologinfo_filtered/annotation_orthologs.tsv"
    shell:
        "python3 scripts/annotation_orthologs_inference.py -f {input} -p {params.percent_identity} -m {params.cdhit_mem} -o {params.output_folder}"