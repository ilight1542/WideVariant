# Annotate assembly using prokka
rule prokka:
    input:
        rules.spades.output.fasta,
    params:
        outdir="results/assembly/prokka/{sampleID}",
    threads: 16
    output:
        txt="results/assembly/prokka/{sampleID}/prokka_out.txt",
        faa="results/assembly/prokka/{sampleID}/prokka_out.faa",
    conda:
        "../envs/prokka.yml"
    shell:
        "prokka --compliant --force --cpus {threads} --outdir {params.outdir} --prefix prokka_out {input}"
