#options for testing:

"""
Additional testing step in snakemake pipeline, then runs tests on all outputs

"""

from ../scripts/etc/bowtie2qc import main

def test_no_outgroup()
    main(input_files_for_this_test)
    assert(files that are output have all 5 samples)
    assert(files that are output match the input)
    assert

def test_all_outgroup():
    main(input_files_for_this_test)

def test_missing_bowtie_file_for_one_sample():
    main(input_files_for_this_test)


input_files_for_one_missing_sample=['my_input_file1.txt','my_input_file3.txt']