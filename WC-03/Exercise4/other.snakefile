from os import path
from os.path import join

configfile: "config.yaml"

workdir: "/students/2020-2021/Thema11/njhindriks/"
FASTQ_DIR = 'samples/'

rule mapping:
    input:
        "genome.fa",
        join(FASTQ_DIR, "{sample}.fastq")
    benchmark:
        "benchmarks/{sample}.benchmark.txt"
    output:
        "mapped_reads/{sample}.bam"
    message: 
        "executing bwa mem on the following {input} to generate the following {output}"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    message:
        "executing samtools sort on the following {input} to generate the following {output}"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    message:
        "executing samtools index on the following {input} to generate the following {output}"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa= "genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    message:
        "executing variant calling with samtools and bcftools with {input}, resulting in {output}" 
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule report:
    input:
        "calls/all.vcf"
    output:
        "out.html"
    message:
        "generating report with {input} resulting in {output}"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeas reference genome 
        and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], metadata="Author: Mr Pipeline", T1=input[0])