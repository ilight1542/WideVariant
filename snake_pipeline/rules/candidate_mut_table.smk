# Builds candidate mutation table (stats across candidate SNV positions)
# Option to build raw coverage matrix and normalized coverage matrix automatically parsed
rule candidate_mutation_table:
    input:
        positions = rules.combine_positions.output.allpositions, # "2-Case/temp/allpositions.pickle",
        string_diversity = rules.candidate_mutation_table_prep.output.string_diversity, # "2-Case/temp/string_diversity_mat.txt",
        string_quals = rules.candidate_mutation_table_prep.output.string_quals, # "2-Case/temp/string_qual_mat.txt",
        string_sampleID_names = rules.candidate_mutation_table_prep.output.string_sampleID_names, # "2-Case/temp/string_sampleID_names.txt",
        string_outgroup_bool = rules.candidate_mutation_table_prep.output.string_outgroup_bool, # "2-Case/temp/string_outgroup_bool.txt",
    output:
        cmt = "2-Case/candidate_mutation_table/group_{cladeID}_candidate_mutation_table.npz",
        cov_norm = "2-Case/candidate_mutation_table/group_{cladeID}_coverage_matrix_norm.npz" if GENERATE_NORMALIZED_COVERAGE_MATRIX == True else [],
        cov_raw = "2-Case/candidate_mutation_table/group_{cladeID}_coverage_matrix_raw.npz" if GENERATE_RAW_COVERAGE_MATRIX == True else [],
    conda:
        "envs/py_for_snakemake.yaml",
    shell:
        "python3 {SCRIPTS_DIRECTORY}/build_candidate_mutation_table.py -p {input.positions} -s {input.string_sampleID_names} -g {input.string_outgroup_bool} -q {input.string_quals} -d {input.string_diversity} -o {output.cmt} -c {output.cov_raw} -n {output.cov_norm} ;"