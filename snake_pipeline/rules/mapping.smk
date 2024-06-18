# Aligns reads to the reference genome with bowtie2
rule bowtie2:
    input:
        fq1 = rules.sickle.output.fq1o,
        fq2 = rules.sickle.output.fq2o,
        bowtie2idx = ancient(rules.refGenome_bowtie_index.output.bowtie2idx) # put here, so rule bowtie2 only executed after rule refGenome_bowtie_index done
    params:
        refGenome = REF_GENOME_DIRECTORY+"/{reference}/genome_bowtie2",
    output:
        samA = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned.sam",
        log = "results/1-mapping/bowtie2/bowtie2_{sampleID}_ref_{reference}.txt", # necessary for bowtie2qc
    conda:
        "../envs/bowtie2.yaml"
    shell:
        # 8 threads coded into json
        "bowtie2 --threads 8 -X 2000 --no-mixed --no-unal --dovetail -x {params.refGenome} -1 {input.fq1} -2 {input.fq2} -S {output.samA} 2> {output.log} "

# Compresses SAM file into BAM file (and removes duplicate reads)
rule sam2bam:
    input:
        samA = rules.bowtie2.output.samA,
    output:
        bamA = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned.sorted.bam",
        bamAidx = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned.sorted.bam.bai",
        DupStats = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_markdup_stats.txt",
    params:
        bamDup = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned_dups.bam",
        bamDupMate = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned_dups.mates.bam",
        bamDupMateSort = "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_aligned_dups.sorted.mates.bam",
    conda:
        "../envs/samtools_bcftools_sm.yaml",
    shadow: # avoids leaving leftover temp files esp if job aborted
        "minimal", 
    shell:
        " samtools view -bS {input.samA} | samtools sort -n - -o {params.bamDup} ;"
        " samtools fixmate -m {params.bamDup} {params.bamDupMate} ;"
        " samtools sort -o {params.bamDupMateSort} {params.bamDupMate} ;"
        " samtools markdup -r -s -f {output.DupStats} -d 100 -m s {params.bamDupMateSort} {output.bamA} ;"
        " samtools index {output.bamA} ;"

# Deletes SAM file once BAM file is created (SAM files are very large)
# must be a to separate rule from sam2bam because rm in sam2bam only deletes link in shadow directory
rule sam2bam_cleanup:
    input:
        bamAidx = rules.sam2bam.output.bamAidx,
    params:
        samA = rules.bowtie2.output.samA,
    output:
        "results/1-mapping/bowtie2/{sampleID}_ref_{reference}_cleanup_done.txt", 
    priority: 100, # prioritizes this rule to get rid of big sam files as fast as possible; default priority for other rules is 0
    shell:
        # -f for cases where sam file doesn't exist (e.g. job previously cancelled/stalled after file deleted but before log file written)
        " rm -f {params.samA} ;"
        " touch {output} ;" 

