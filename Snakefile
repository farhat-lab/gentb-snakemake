#configfile: "config.yaml"
(SAMPLES,) = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

#reference = "test-data/reference/NC_000962.3_Mycobacterium tuberculosis_H37Rv.fa"

#or 
#wildcards = glob_wildcards("fastq/{sample}_R1.fastq.gz")
#then use {wildcards.sample}


rule all:
    input:
        expand("results/{sample}.quality-control.html", sample=SAMPLES)

rule fastqc:
    input: 
        fwd="data/fastq/{sample}_R1.fastq.gz", rev="data/fastq/{sample}_R2.fastq.gz"    
    output:
        html="results/{sample}.quality-control.html",
        zip="results/{sample}_quality-control_fastqc.zip" # in case we want to use multiqc later as it requires this suffix
    log:
        "logs/fastqc/{sample}.log"
    message: 
        "validating FASTQ files"
    wrapper:
        "0.36.0/bio/fastqc"


# rule fastptrim:
#     input:
#         fwd="test-data/{sample}_R1.fastq.gz", rev="test-data/{sample}_R2.fastq.gz"
#     output:
#         fwd="results/{sample}_R1_trimmed.fastq", rev="results/{sample}_R2_trimmed.fastq "
#         htmlout="results/{sample}_fastp.html", jsonout="results/{sample}_fastp.json"
#     message: 
#         "* trimming reads (fastp)"
#     shell: 
#         "fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} --html={output.htmlout} --json={output.jsonout}"


# #Kraken2 Classify
# module load gcc/6.2.0 module load python/3.6.0 kraken2_wrapper.sh 
# ${file}_trimmed_1.fastq ${file}_trimmed_2.fastq @{file}_kraken2_classified.txt 
# @{file}_kraken2_report.txt