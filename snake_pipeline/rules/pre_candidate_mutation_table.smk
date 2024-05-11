# Generates a list of candidate SNV positions for a given sample
rule variants2positions:
    input:
        variants = "results/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz",
    params:
        refGenomeDir = REF_GENOME_DIRECTORY+"/{reference}/",
        outgroup_tag = "{outgroup}", # boolean (0==ingroup or 1==outgroup)
        maxFQ = -30,
    output:
        positions = "results/2-case/temp/{sampleID}_ref_{reference}_outgroup{outgroup}_positions.npz",
    conda:
        "../envs/py_for_snakemake.yaml",
    shell:
        "mkdir -p results/2-case/temp/ ;"
        "python {SCRIPTS_DIRECTORY}/variants2positions.py -i {input.variants} -o {output.positions} -r {params.refGenomeDir} -q {params.maxFQ} -b {params.outgroup_tag} ;"    


# Creates a list of files with candidate SNV positions from each sample
rule combine_positions_prep:
    input:
        positions = get_positions_prep,
    output:
        string_input_p_positions = "results/2-case/temp/group_{cladeID}_string_file_other_p_to_consider.txt",
    run:
        with open( output.string_input_p_positions ,"w") as f: 
            print(*input.positions, sep="\n", file=f)


# Build input for candidate_mutation_table
rule candidate_mutation_table_prep:
    input:
        diversity = get_diversity,
        quals = get_quals,
    params:
        sampleID_names = get_sampleID_names,
        outgroup_bool = get_outgroup_bool,
    output:
        string_diversity = "results/2-case/temp/group_{cladeID}_string_diversity.txt",
        string_quals = "results/2-case/temp/group_{cladeID}_string_qual.txt",
        string_sampleID_names = "results/2-case/temp/group_{cladeID}_string_sampleID_names.txt",
        string_outgroup_bool = "results/2-case/temp/group_{cladeID}_string_outgroup_bool.txt",
    run:
        with open( output.string_diversity ,"w") as f: 
            print(*input.diversity, sep="\n", file=f)
        with open( output.string_quals ,"w") as f:
            print(*input.quals, sep="\n", file=f)
        with open( output.string_sampleID_names ,"w") as f: 
            print(*params.sampleID_names, sep="\n", file=f)
        with open( output.string_outgroup_bool ,"w") as f: 
            print(*params.outgroup_bool, sep="\n", file=f)


# Generates a list of candidate SNV positions based on candidate SNV positions across ingroup samples
# (Ignores candidate SNVs in samples marked as outgroups)
rule combine_positions:
    input:
        string_input_pos = rules.combine_positions_prep.output.string_input_p_positions,
        string_outgroup_bool = rules.candidate_mutation_table_prep.output.string_outgroup_bool,
    params:
        # file_other_p_to_consider = "results/2-case/temp/other_positions.npz",
        refGenomeDir = get_ref_genome, # expands to single reference genome!
    output:
        allpositions = "results/2-case/temp/group_{cladeID}_allpositions.npz",
    conda:
        "../envs/py_for_snakemake.yaml",
    shell:
        "python {SCRIPTS_DIRECTORY}/combine_positions.py -i {input.string_input_pos} -r {params.refGenomeDir} -b {input.string_outgroup_bool} -o {output.allpositions} ;"