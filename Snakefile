configfile: "config.yaml"
reference=config["reference"]

(SAMPLES,) = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("results/{sample}_R1_trimmed.fastq", sample=SAMPLES),
        expand("results/{sample}_aln.sam", sample=SAMPLES),
        expand("results/{sample}.namesorted.bam", sample=SAMPLES),
        expand("results/{sample}.fixmate.bam", sample=SAMPLES),
        expand("results/{sample}.positionsort.bam", sample=SAMPLES),
        expand("results/{sample}.rmdup.bam", sample=SAMPLES),
        expand("results/{sample}.rmdup.bam.bai", sample=SAMPLES),
        expand("results/{sample}.vcf", sample=SAMPLES)

rule fastp:
    input:
        sample=["data/fastq/{sample}_R1.fastq.gz", "data/fastq/{sample}_R2.fastq.gz"]
    output:
        trimmed=["results/{sample}_R1_trimmed.fastq", "results/{sample}_R2_trimmed.fastq"],
        html="results/{sample}.fastp.html",
        json="results/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        extra=""
    threads: 2
    wrapper:
        "0.36.0/bio/fastp"

#cite 10.1093/bioinformatics/bty560

#insert here Kraken2 Classify

rule minimap2:
    input:
        fwd="results/{sample}_R1_trimmed.fastq", rev="results/{sample}_R2_trimmed.fastq"
    output: 
        temp("results/{sample}_aln.sam")
    threads: 3
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 -ax sr {reference} {input.fwd} {input.rev} > {output}"

#cite 10.1093/bioinformatics/bty191

# use alternatively the minimap2 wrapper as below
# rule minimap2:
#     input:
#         target, 
#         query=["data/results/{sample}_R1_trimmed.fastq", "data/results/{sample}_R2_trimmed.fastq"]
#     output:
#         "results/{sample}_aln.sam"
#     log:
#         "logs/minimap2/{sample}.log"
#     params:
#         extra="-ax sr"
#     threads: 3
#     wrapper:
#         "0.36.0/bio/minimap2/aligner"

rule samtools_sort:
    input: 
        "results/{sample}_aln.sam"
    output: 
        "results/{sample}.namesorted.bam"
    shell: 
        "samtools sort --threads 8 -n -m 4G -o {output} {input}"

rule samtools_fixmate:
    input: 
        "results/{sample}.namesorted.bam"
    output:
        "results/{sample}.fixmate.bam"
    shell:
        "samtools fixmate -m {input} {output}"

rule samtools_positionsort:
    input:
        "results/{sample}.fixmate.bam"
    output:
        "results/{sample}.positionsort.bam"
    shell:
        "samtools sort -o {output} {input}"

rule samtools_rmdup:
    input: 
        "results/{sample}.positionsort.bam"
    output:
        "results/{sample}.rmdup.bam"
    shell:
        "samtools markdup --threads 4 -r -s {input} {output}"

#insert here "check depth" using the python script checkd.py

rule samtools_index:
    input:
        "results/{sample}.rmdup.bam"
    output:
        "results/{sample}.rmdup.bam.bai"
    shell:
        "samtools index {input} {output}"


#cite 10.1093/bioinformatics/btp352

rule pilon:
    input: 
        "results/{sample}.rmdup.bam"
    output: 
        "results/{sample}.vcf"
    shell:
        "java -Xmx12G -jar pilon --genome {reference} --bam {input} --output {output} --vcf"

#cite 10.1371/journal.pone.0112963

#check here on how Pilon is called: /Users/matthias/miniconda3/envs/snakemake-gentb/bin/pilon. Specify the amount of RAM is used.
#see lucas megapipe for it



