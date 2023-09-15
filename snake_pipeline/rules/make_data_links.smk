rule make_data_links:
    # NOTE: All raw data needs to be names fastq.gz. No fq! The links will be names fq though.
    input:
        sample_info_csv = "data/{sampleID}/sample_info.csv",
    output:
        # Recommend using symbolic links to your likely many different input files
        fq1 = "data/{sampleID}/R1.fq.gz",
        fq2 = "data/{sampleID}/R2.fq.gz",
    run:
    # create links
        paths, sample, reference, filename = read_sample_info_CSV(input.sample_info_csv)
        if len(paths) > 1:
            cp_append_files(paths,sample,filename,'data')
        else:
            makelink(paths[0],sample,filename,'data')