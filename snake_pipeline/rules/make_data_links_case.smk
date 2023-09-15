rule make_data_links_case:
    # input:
    #     sample_info_csv="data/{sampleID}/sample_info.csv",
    output:
        # Recommend using symbolic links to your likely many different input files      
        vcf_links = expand("1-Mapping/vcf/{sampleID}_ref_{references}_aligned.sorted.strain.variant.vcf.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls),
        qual_links = expand("1-Mapping/quals/{sampleID}_ref_{references}_outgroup{outgroup}.quals.pickle.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
        div_links = expand("1-Mapping/diversity/{sampleID}_ref_{references}_outgroup{outgroup}.diversity.pickle.gz",zip,sampleID=SAMPLE_ls, references=REF_Genome_ls, outgroup=OUTGROUP_ls),
    run:
        subprocess.run( "mkdir -p 1-Mapping/vcf/ 1-Mapping/quals/ 1-Mapping/diversity/ " ,shell=True)
        for idx, ele in enumerate(SAMPLE_ls):
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/1-Mapping/diversity/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*diversity* 1-Mapping/diversity/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_outgroup{OUTGROUP_ls[idx]}.diversity.pickle.gz" ,shell=True)
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/1-Mapping/quals/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*quals* 1-Mapping/quals/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_outgroup{OUTGROUP_ls[idx]}.quals.pickle.gz" ,shell=True)
            subprocess.run( f"ln -fs -T {PATH_ls[idx]}/1-Mapping/vcf/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_*variant.vcf.gz 1-Mapping/vcf/{SAMPLE_ls[idx]}_ref_{REF_Genome_ls[idx]}_aligned.sorted.strain.variant.vcf.gz" ,shell=True)