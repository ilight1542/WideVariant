rule make_data_links_case:
    # input:
    #     sample_info_csv="results/data/{sampleID}/sample_info.csv",
    output:
        # Recommend using symbolic links to your likely many different input files      
        vcf_links = expand("results/1-mapping/vcf/{sampleID}_ref_{references}_aligned.sorted.strain.variant.vcf.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls),
        qual_links = expand("results/1-mapping/quals/{sampleID}_ref_{references}_outgroup{outgroup}.quals.pickle.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
        div_links = expand("results/1-mapping/diversity/{sampleID}_ref_{references}_outgroup{outgroup}.diversity.pickle.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
    run:
        subprocess.run( "mkdir -p results/1-mapping/vcf/ results/1-mapping/quals/ results/1-mapping/diversity/ " ,shell=True)
        for idx, ele in enumerate(SAMPLE_ls):
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/results/1-mapping/diversity/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*diversity* results/1-mapping/diversity/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_outgroup{OUTGROUP_ls[idx]}.diversity.pickle.gz" ,shell=True)
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/results/1-mapping/quals/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*quals* results/1-mapping/quals/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_outgroup{OUTGROUP_ls[idx]}.quals.pickle.gz" ,shell=True)
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/results/1-mapping/vcf/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*variant.vcf.gz results/1-mapping/vcf/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_aligned.sorted.strain.variant.vcf.gz" ,shell=True)
