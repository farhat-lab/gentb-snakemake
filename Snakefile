configfile: "data/config.yaml"
reference=config["reference"]
reference_fasta=config["reference_fasta"]
lineage=config["lineage"]
genomesummary=config["genomesummary"]
noncodingsummary=config["noncodingsummary"]
variant_name_list=config["variant_name_list"]
variant_name_list_pyrazinamide=config["variant_name_list_pyrazinamide"]
lineage_snp_mf=config["lineage_snp_mf"]

(SAMPLES,)=glob_wildcards("data/fastq/{sample}_1.fastq.gz")

ruleorder: samtools_index > pilon > annotator > generate_matrix > TBpredict_all > enhance

onstart:
    "Rscript scripts/initialise.R"

rule all:
    input:
        expand("results/{sample}.final.predict.json", sample=SAMPLES)

rule fastp:
    input:
        fwd="data/fastq/{sample}_1.fastq.gz", 
        rev="data/fastq/{sample}_2.fastq.gz",
    output:
        fwd=temp("results/{sample}_1_trimmed.fastq.gz"), 
        rev=temp("results/{sample}_2_trimmed.fastq.gz"),
    shell:
        "fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev}"

rule minimap2:
    input:
        fwd="results/{sample}_1_trimmed.fastq.gz", 
        rev="results/{sample}_2_trimmed.fastq.gz"
    output: 
        "results/{sample}_aln.sam"
    shell:
        "minimap2 -t 1 -ax sr {reference} {input.fwd} {input.rev} > {output}"

#during installation make sure to use `conda install -c bioconda samtools=1.9`
rule samtools_sort:
    input: 
        "results/{sample}_aln.sam"
    output: 
        "results/{sample}.sorted.bam"
    shell: 
        "samtools sort -m 15G -n -o {output} {input}"


rule samtools_fixmate:
    input:
        "results/{sample}.sorted.bam"
    output:
        temp("results/{sample}.fixmate.bam")
    shell:
        "samtools fixmate -m {input} {output}"

rule samtools_positionsort:
    input:
        "results/{sample}.fixmate.bam"
    output:
        temp("results/{sample}.positionsort.bam")
    shell:
        "samtools sort -o {output} {input}"

rule samtools_rmdup:
    input: 
        "results/{sample}.positionsort.bam"
    output:
        "results/{sample}.rmdup.bam"
    shell:
        "samtools markdup -r -s {input} {output}"


rule samtools_index:
    input:
        "results/{sample}.rmdup.bam"
    output:
        "results/{sample}.rmdup.bam.bai"
    shell:
        "samtools index -b {input} {output}"

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
        "python3 scripts/vcf_cutter.py -i {input} -o {output}"

rule annotator:
    input:
        "results/{sample}.cut.vcf"
    output:
        "results/{sample}.var"
    shell:
        "perl scripts/flatAnnotatorVAR_2.0.1.pl {reference_fasta} {genomesummary} {noncodingsummary} {input} 10 0.4 PASS AMB > {output}"

rule generate_matrix:
    input:
        "results/{sample}.var"
    output:
        "results/{sample}.matrix.csv"
    shell:
        "python3 scripts/generate_matrix.py {variant_name_list} {input} > {output}"

rule generate_matrix_pyrazinamide:
    input:
        "results/{sample}.var"
    output:
        "results/{sample}.matrix.pyrazinamide.csv"
    shell:
        "python3 scripts/generate_matrix_pyrazinamide.py {variant_name_list_pyrazinamide} {input} > {output}"
        
rule TBpredict_all:
    input:
        RF="results/{sample}.matrix.csv", PZA="results/{sample}.matrix.pyrazinamide.csv"
    output: 
        RF="results/{sample}.matrix.json", PZA="results/{sample}.matrix.pza.json"
    shell:
        "Rscript scripts/TBpredict_combined.R {input.RF} {input.PZA} && ls {output.RF} {output.PZA}"
 
rule merge_pyrazinamide_prediction:
    input:
        all_predictions="results/{sample}.matrix.json",
        pza_prediction="results/{sample}.matrix.pza.json",
    output:
        "results/{sample}.predict.json"
    shell:
        "python3 scripts/merge_retrained_pyrazinamide_prediction.py {input.all_predictions} {input.pza_prediction} {output}" 

rule enhance:
    input: 
        matrix="results/{sample}.predict.json", var="results/{sample}.var"
    output: 
        "results/{sample}.final.predict.json"
    shell: 
        "python3 scripts/varMatchUnk.py {input.var} {lineage_snp_mf} {input.matrix} -o {output}"


