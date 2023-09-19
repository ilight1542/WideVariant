# Processes BAM file into VCF files
rule mpileup2vcf:
    input:
        bamA = rules.sam2bam.output.bamA,
        bamClean = rules.sam2bam_cleanup.output,
        fasta_idx = ancient(rules.samtools_idx.output.fasta_idx),
    params:
        ref = REF_GENOME_DIRECTORY+"/{reference}/genome.fasta",
        vcf_raw = "results/results/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.gz",
    output:
        pileup = "results/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.pileup",
        variants = "results/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.variant.vcf.gz",
        vcf_strain = "results/1-mapping/vcf/{sampleID}_ref_{reference}_aligned.sorted.strain.vcf.gz",
    conda:
        "../envs/samtools_bcftools_sm.yaml"
    shadow: 
        "shallow", # avoids leaving leftover temp files esp if job aborted, but need shallow to allow access to refgenome folder
    shell:
        " samtools mpileup -q30 -x -s -O -d3000 -f {params.ref} {input.bamA} > {output.pileup} ;" 
        " samtools mpileup -q30 -t SP -d3000 -vf {params.ref} {input.bamA} > {params.vcf_raw} ;"
        " bcftools call -c -Oz -o {output.vcf_strain} {params.vcf_raw} ;"
        " bcftools view -Oz -v snps -q .75 {output.vcf_strain} > {output.variants} ;"
        " tabix -p vcf {output.variants} ;"
        " rm {params.vcf_raw}"

# Parses VCF with python script
rule vcf2quals:
    input:
        vcf_strain = rules.mpileup2vcf.output.vcf_strain,
    params:
        refGenomeDir = REF_GENOME_DIRECTORY+"/{reference}/",
    output:
        file_quals = "results/1-mapping/quals/{sampleID}_ref_{reference}_outgroup{outgroup}.quals.npz",
    conda:
        "../envs/py_for_snakemake.yaml",
    shell:
        "mkdir -p results/1-mapping/quals/ ;"
        "python {SCRIPTS_DIRECTORY}/vcf2quals_snakemake.py -i {input.vcf_strain} -r {params.refGenomeDir} -o {output.file_quals} ;"

# Parses pileup with python script
rule pileup2diversity_matrix:
    input:
        pileup = rules.mpileup2vcf.output.pileup,
    params:
        refGenomeDir = REF_GENOME_DIRECTORY+"/{reference}/", 
    output:
        file_diversity = "results/1-mapping/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.diversity.npz",
        file_coverage = "results/1-mapping/diversity/{sampleID}_ref_{reference}_outgroup{outgroup}.aligned.sorted.strain.variant.coverage.npz",
    conda:
        "../envs/py_for_snakemake.yaml",
    shell:
        "mkdir -p results/1-mapping/diversity/ ;"
        "python {SCRIPTS_DIRECTORY}/pileup2diversity.py -i {input.pileup} -r {params.refGenomeDir} -o {output.file_diversity} -c {output.file_coverage} ;"