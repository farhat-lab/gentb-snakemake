## GenTB-RandomForest tuberculosis resistance prediction pipeline in Snakemake

This repository contains the analysis pipeline used by [GenTB](https://gentb.hms.harvard.edu) to predict antibiotic resistance _M. tuberculosis_ wrapped in _[Snakemake](https://snakemake.readthedocs.io/en/stable/)_. The code for the [GenTB webiste](https://gentb.hms.harvard.edu) is described [here](https://github.com/farhat-lab/gentb-site).

### Overview

The pipeline uses well established programs to align reads to a reference strain and call variants. Raw sequence reads, in the form of a fastq file, are aligned to H37Rv reference strain using minimap2. Samtools is employed for sorting of the alignment file and removal of duplicate reads. Variants are called with Pilon. The variants are then annotated and filtered using custom scripts and fed to two multivariate prediction models ([GenTB-Random Forest](https://www.atsjournals.org/doi/full/10.1164/rccm.201510-2091OC) and [GenTB Wide and Deep Neural Network](https://www.sciencedirect.com/science/article/pii/S2352396419302506?via%3Dihub)). 

<img src="https://github.com/mgro/mgro.github.io/blob/master/GenTB%20Pipeline.png" width="700" height="700" />



### Note

This `Snakemake` implementation of GenTB-Random Forest currently does not include the Kraken and Check Depth steps.

### Usage

1. Clone this Github repository to your local machine

```
git clone https://github.com/farhat-lab/gentb-snakemake.git
```

2. Download example fastq files:

```
#change directory into the repository folder, create a directory to save the fastq files, then download sample fastq files
cd gentb-snakemake/; mkdir data/fastq; cd data/fastq/

#download one set of paired-ended fastQ files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/009/ERR1873539/ERR1873539_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/009/ERR1873539/ERR1873539_2.fastq.gz

#change back to main `gentb-snakemake` directory
cd ../../
```

2. Make sure you have the package manager [conda](https://docs.conda.io/en/latest/miniconda.html) installed to source all software and dependencies into a new environment:
```
conda env create -f envs/snakemake-gentb.yml
conda activate snakemake-gentb
```

4. Start the analysis pipeline in _Snakemake_
```
snakemake
```

5. Start the analysis pipeline in _Snakemake_ on a SLURM cluster:
```
snakemake -s Snakefile -j 100 --cluster-config cluster.json --cluster "sbatch --mem {cluster.mem} -t {cluster.t} -c {cluster.c} -p {cluster.p}"
```

### Output description

The final output file `<sampleID>.predict.json` is printed in JSON format and can be parsed using R or python functions. The output consists of three parts.

1. Drug resistance probabilities

A list of drug-lists where each drug-list has the following items:
```
["<sampleID>", 
 "<drug>", 
 "<binary resistance (1/0)>", 
 "<False Negative Rate>", 
 "<False Positive Rate>", 
 "<probability of resistance>"]
 ```
 
 For the example of `ERR1873539` and the drug ethambutol, this part of the output would look as follows:
 ```
 ["ERR1873539", 
  "emb", 
  "1", 
  "21.28", 
  "16.64", 
  "0.998"]
  ```
  
  2. Important drug resistance mutations
  
Next, a list of five lists of important drug resistance mutations (defined [here](https://www.atsjournals.org/doi/full/10.1164/rccm.201510-2091OC)). The items per list are ordered the same like the resistance probability lists above. 
  
For the example of `ERR1873539`, this part of the output would look as follows:
  
  ```
"ERR1873539": [
[
  "SNP_CN_2155168_C944G_S315T_katG", 
  "SNP_CN_761110_A1304T_D435V_rpoB", 
  "INS_CF_2288725_i517G_173G_pncA", 
  "SNP_CN_4247431_G918A_M306I_embB", 
  "SNP_N_1472359_A514C_rrs", 
  "SNP_CN_4326333_C1141G_A381P_ethA", 
  "SNP_N_1473246_A1401G_rrs", 
  "SNP_N_1473246_A1401G_rrs", 
  "SNP_N_1473246_A1401G_rrs", 
  "SNP_CN_7582_A281G_D94G_gyrA", 
  "SNP_CN_7582_A281G_D94G_gyrA", 
  "SNP_CN_7582_A281G_D94G_gyrA", 
  null
], 
[
  "SNP_CN_4247431_G918A_M306I_embB", 
  null, 
  null, 
  null, 
  "SNP_N_1473246_A1401G_rrs", 
  "SNP_P_1673423_G17T_promoter_fabG1.inhA", 
  null, 
  null, 
  "SNP_N_1472359_A514C_rrs", 
  null, 
  null, 
  null, 
  null
], 
[
  "SNP_P_1673423_G17T_promoter_fabG1.inhA", 
  null, 
  null, 
  null, 
  "SNP_CN_4407927_T276G_E92D_gid", 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null
], 
[
  null, 
  null, 
  null, 
  null, 
  "SNP_CN_4407967_A236G_L79S_gid", 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null
], 
[
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null, 
  null
 ]
]


3. Last, another list of five list follows with other or novel variants that were not seen durinig training the Random Forest model but occurred in known resistance genes. 


Tested on CentOS Linux Version 7
