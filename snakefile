configfile: "config.yaml"
reference=config["reference"]
reference_fasta=config["reference_fasta"]
lineage=config["lineage"]
genomesummary=config["genomesummary"]
noncodingsummary=config["noncodingsummary"]
variant_name_list=config["variant_name_list"]
lineage_snp_mf=config["lineage_snp_mf"]

(SAMPLES,) = glob_wildcards("/n/scratch2/mg451/fastq/all-fastq-01/{sample}_1.fastq.gz")

ruleorder: samtools_index > pilon > annotator > generate_matrix > TBpredict > enhance

onstart:
    "Rscript ./scripts/initialise.R"

rule all:
    input:
        #expand("/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_1_trimmed.fastq.gz", sample=SAMPLES),
        #expand("/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_aln.sam", sample=SAMPLES),
        #expand("results/{sample}.sorted.bam", sample=SAMPLES),
        #expand("results/{sample}.rmdup.bam", sample=SAMPLES),
        #expand("results/{sample}.check_depth_fake_file", sample=SAMPLES),
        #expand("results/{sample}.rmdup.bam.bai", sample=SAMPLES),
        #expand("results/{sample}.vcf", sample=SAMPLES),
        #expand("results/{sample}.cut.vcf", sample=SAMPLES),
        #expand("results/{sample}.lineage.txt", sample=SAMPLES),
        expand("results/{sample}.var", sample=SAMPLES),
        expand("results/{sample}.matrix.csv", sample=SAMPLES), 
        expand("results/{sample}.matrix.json", sample=SAMPLES),
        expand("results/{sample}.final.matrix.json", sample=SAMPLES)

rule fastp:
    input:
        fwd="/n/scratch2/mg451/fastq/all-fastq-01/{sample}_1.fastq.gz", 
        rev="/n/scratch2/mg451/fastq/all-fastq-01/{sample}_2.fastq.gz"
    output:
        fwd=temp("/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_1_trimmed.fastq.gz"), 
        rev=temp("/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_2_trimmed.fastq.gz")
    log:
        "logs/fastp/{sample}.log"
    params:
        extra=""
    threads: 1
    shell:
        "fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev}"

#cite 10.1093/bioinformatics/bty560

rule minimap2:
    input:
        fwd="/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_1_trimmed.fastq.gz", 
        rev="/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_2_trimmed.fastq.gz"
    output: 
        temp("/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_aln.sam")
    threads: 1
    log:
        "logs/minimap2/{sample}.log"
    shell:
        "minimap2 -t 1 -a {reference} {input.fwd} {input.rev} > {output}"

#cite 10.1093/bioinformatics/bty191

rule samtools_sort:
    input: 
        "/n/scratch2/mg451/gentb_benchmarking_temp/{sample}_aln.sam"
    output: 
        temp("results/{sample}.sorted.bam")
    shell: 
        "samtools sort -m 15G -o {output} {input}"

rule samtools_rmdup:
    input: 
        "results/{sample}.sorted.bam"
    output:
        "results/{sample}.rmdup.bam"
    shell:
        "samtools markdup -r -s {input} {output}"

rule depthcheck:
    input: 
        "results/{sample}.rmdup.bam"
    output:
        "results/{sample}.check_depth_fake_file"
    shell:
        "samtools depth -a {input} | python3 ./scripts/checkd.py; "
        "touch {output}"


rule samtools_index:
    input:
        "results/{sample}.rmdup.bam"
    output:
        "results/{sample}.rmdup.bam.bai"
    shell:
        "samtools index -b {input} {output}"

#cite 10.1093/bioinformatics/btp352

rule pilon:
    input: 
        bam="results/{sample}.rmdup.bam", bai="results/{sample}.rmdup.bam.bai"
    output: 
        temp("results/{sample}.vcf"),
        temp("results/{sample}.fasta")
    params: 
        jar_path=lambda wildcards, output: os.path.join(os.environ["CONDA_PREFIX"], 'share', 'pilon-1.23-2'),
        output=lambda wildcards, output: os.path.join(output[0].rsplit('/', 1)[0], wildcards[0])
    shell:
        "java -Xmx15G -jar {params.jar_path}/pilon-1.23.jar --threads 1 --variant --genome {reference_fasta} "
        "--bam {input.bam} --output {params.output} --vcf"

rule cut_vcf:
    input:
        "results/{sample}.vcf"
    output:
        "results/{sample}.cut.vcf"
    shell:
        "python3 ./scripts/vcf_cutter.py -i {input} -o {output}"
    
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
        "results/{sample}.cut.vcf"
    output:
        "results/{sample}.var"
    shell:
        "perl ./scripts/flatAnnotatorVAR_2.0.1.pl {reference_fasta} {genomesummary} {noncodingsummary} {input} 10 0.4 PASS AMB > {output}"

rule generate_matrix:
    input:
        "results/{sample}.var"
    output:
        "results/{sample}.matrix.csv"
    shell:
        "python3 ./scripts/generate_matrix.py {variant_name_list} {input} > {output}"

rule TBpredict:
    input:
        "results/{sample}.matrix.csv"
    output: 
        "results/{sample}.matrix.json"
    shell:
        "Rscript ./scripts/TBpredict.R {input}"

rule enhance:
    input: 
        matrix="results/{sample}.matrix.json", var="results/{sample}.var"
    output: 
        "results/{sample}.final.matrix.json"
    shell: 
        "python3 ./scripts/varMatchUnk.py {input.var} {lineage_snp_mf} {input.matrix} -o {output}"


