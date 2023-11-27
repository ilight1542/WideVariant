# Get two-column (caldeID,path2faa) input file for ortholog_inference script
rule build_annotation_orthologs_input:
    input:
        bakta_faa=expand("results/assembly/bakta/{sampleID}/bakta_out.faa",sampleID=SAMPLE_ls),
    params:
        clade_identifier=expand("{sampleID}",sampleID=SAMPLE_ls),
    output:
        "results/assembly/orthologinfo_filtered/input_files.tsv",
    group:
        "infer_orthologs"
    shell:
        """
        paste <(echo {params.clade_identifier} | {SCRIPTS_DIRECTORY}/sed_nl.sh ) <(echo {input.bakta_faa} | {SCRIPTS_DIRECTORY}/sed_nl.sh ) > {output}
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


# Build symbolic links with unique basenames for roary
rule roary_input: 
    input:
        bakta_gff="results/assembly/bakta/{sampleID}/bakta_out.gff3",
    output:
        roary_infile="results/assembly/roary/tmp_input_data/{sampleID}.gff3",
    params:
        startline=9,
    shell:
        " tail -n +{params.startline} $(pwd)/{input.bakta_gff} > {output.roary_infile} " 
    
rule roary:
    input:
        bakta_datalinks=expand("results/assembly/roary/tmp_input_data/{sampleID}.gff3",sampleID=SAMPLE_ls),
    output:
        gpa_csv="results/assembly/roary/results/gene_presence_absence.csv",
    params:
        outdir="results/assembly/roary/results",
        blastp_ident=95,
        mcl_iv=1.5, 
    threads: 4,
    conda: 
        "../envs/roary.yaml",
    shell:
        " rm -rf {params.outdir} ; "
        " roary -p {threads} -f {params.outdir} -v -i {params.blastp_ident} -iv {params.mcl_iv} {input.bakta_datalinks} "

rule filter_roary:
    input:
        gpa_csv=rules.roary.output.gpa_csv,
    output:
        orthologs="results/assembly/roary/annotation_orthologs.tsv",
    params:
        roary_input="results/assembly/roary/tmp_input_data/"
    shell:
        """ rm -rf {params.roary_input} ;"""
        """ awk -F'","' 'BEGIN{{OFS="\t"}} NR==1 {{printf "ortholog_id"; for(i=15; i<=NF; i++){{ gsub("\"", "", $i); printf "\t%s", $i }} printf "\\n"}}' {input.gpa_csv} > {output.orthologs} ;"""
        """ awk -F'","' 'BEGIN{{OFS="\t"; count=0}} NR>1 {{count++; gsub("\"", "", $1); printf "RGOI_%05d", count; for(i=15; i<=NF; i++){{ gsub("\"", "", $i); printf "\t%s", $i }} printf "\\n"}}' {input.gpa_csv} >> {output.orthologs} ;"""
        
 
##TODO: Docs, need to include conda env!!!
