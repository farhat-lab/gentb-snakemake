configfile: "config.yaml"
reference=config["reference"]
lineage=config["lineage"]
genomesummary=config["genomesummary"]
noncodingsummary=config["noncodingsummary"]

(SAMPLES,) = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("results/{sample}_R1_trimmed.fastq", sample=SAMPLES),
        expand("results/{sample}_aln.sam", sample=SAMPLES),
        expand("results/{sample}.namesorted.bam", sample=SAMPLES),
        expand("results/{sample}.fixmate.bam", sample=SAMPLES),
        expand("results/{sample}.positionsort.bam", sample=SAMPLES),
        expand("results/{sample}.rmdup.bam", sample=SAMPLES),
        expand("results/fake_file.{sample}.rmdup.bam", sample=SAMPLES),
        expand("results/{sample}.rmdup.bam.bai", sample=SAMPLES),
        expand("results/{sample}.vcf", sample=SAMPLES),
        expand("results/{sample}.lineage.txt", sample=SAMPLES),
        expand("results/{sample}.var", sample=SAMPLES)

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

# rule kraken2:
#     input:
#     output:

# from Luca
# rule kraken2:
#     input: fwd="results/{sample}_R1_trimmed.fastq", rev="results/{sample}_R2_trimmed.fastq"
#     output: classified="results/{sample}_kraken2_classified.txt", report="results/{sample}_kraken2_report.txt"
#     shell: "kraken2 --threads 4 --db krakenDB/tbdb --paired {input.fwd} {input.fwd} > {output.classified} 2> {output.report} && cat {output.report} 1>&2 && cat {output.report} | python3 bin/classify.py"
#cite 10.1186/gb-2014-15-3-r46

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

rule depthcheck:
    input: 
        "results/{sample}.rmdup.bam"
    output:
        "results/fake_file.{sample}.rmdup.bam"
    shell:
        "samtools depth -a {input} | python3 ./scripts/checkd.py; "
        "touch {output}"

# --keepgoing
# -k

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
        "results/{sample}.vcf",
        temp("results/{sample}.fasta")
    params: 
        jar_path=lambda wildcards, output: os.path.join(os.environ["CONDA_PREFIX"], 'share', 'pilon-1.23-1'),
        output=lambda wildcards, output: os.path.join(output[0].rsplit('/', 1)[0], wildcards[0])
    shell:
        "java -Xmx12G -jar {params.jar_path}/pilon-1.23.jar --genome {reference} "
        "--bam {input} --output {params.output} --vcf"
    
#cite 10.1371/journal.pone.0112963

rule call_lineage:
    input: 
        "results/{sample}.vcf"
    output:
        "results/{sample}.lineage.txt"
    shell:
        "python3 ./scripts/fast-lineage-caller-vcf.py {input} {lineage} > {output}"

rule annotator:
    input:
        "results/{sample}.vcf"
    output:
        "results/{sample}.var"
    shell:
        "perl ./scripts/flatAnnotatorVAR_2.0.pl {reference} {genomesummary} {noncodingsummary} {input} 10 0.4 PASS AMB > {output}"
