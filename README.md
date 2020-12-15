## GenTB FASTQ analysis pipeline

This repository contains the analysis pipeline used by [GenTB](https://gentb.hms.harvard.edu) to predict antibiotic resistance _M. tuberculosis_ wrapped in _[Snakemake](https://snakemake.readthedocs.io/en/stable/)_. The code for the [GenTB webiste](https://gentb.hms.harvard.edu) is described [here](https://github.com/farhat-lab/gentb-site).

### Overview

The pipeline uses well established programs to align reads to a reference strain and call variants. Raw sequence reads, in the form of a fastq file, are aligned to H37Rv reference strain using minimap2. Samtools is employed for sorting of the alignment file and removal of duplicate reads. Variants are called with Pilon. The variants are then annotated and filtered using custom scripts and fed to two multivariate prediction models ([GenTB-Random Forest](https://github.com/mgro/mgro.github.io/blob/master/GenTB%20Pipeline.png) and [GenTB Wide and Deep Neural Network](https://www.sciencedirect.com/science/article/pii/S2352396419302506?via%3Dihub)). 

<img src="https://github.com/mgro/mgro.github.io/blob/master/GenTB%20Pipeline.png" width="700" height="700" />



### Usage

1. Clone this Github repository to your local machine

2. Make sure you have the package manager [conda](https://docs.conda.io/en/latest/miniconda.html) installed to source all software and dependencies into a new environment:
```
conda env create -f envs/snakemake-gentb.yml
conda activate snakemake-gentb
```

3. Download example fastq files:

```
mkdir data/fastq; cd data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/009/ERR1873539/ERR1873539_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/009/ERR1873539/ERR1873539_2.fastq.gz
```

4. Start the analysis pipeline in _Snakemake_
```
snakemake -s snakefile_test_run --use-conda
```

5. Start the analysis pipeline in _Snakemake_ on a SLURM cluster:
```
snakemake -s snakefile_test_run -j 100 --cluster-config cluster.json --cluster "sbatch --mem {cluster.mem} -t {cluster.t} -c {cluster.c} -p {cluster.p}"
```
