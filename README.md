# nf-gwas-pipeline
A Nextflow Genome-Wide Association Study (GWAS) Pipeline

### Clone Repository
```bash

$ git clone https://github.com/montilab/nf-gwas-pipeline

```
### Initalize Paths to Test Data

We have provided multiple toy datasets for testing the pipeline and ensuring all paths and dependencies are properly setup. To set the toy data paths to your local directory, run the following script.

```bash
$ cd nf-gwas-pipeline
$ python utils/paths.py  
```

### Download Nextflow Executable
```
$ curl -s https://get.nextflow.io | bash
```

### Quick Start with Docker

We have created a pre-built Docker image with all of the dependencies installed. To get started, first make sure [Docker is installed](https://docs.docker.com/get-docker/). Then pull down the image onto your local machine.

```
$ docker pull montilab/gwas:latest
```

Optionally you could build this image  yourself from the Dockerfile which specifies all of the dependencies required. *Note: This might take a while!*

```
$ docker build --tag montilab/gwas:latest .
```

### Run with docker
```
$ ./nextflow gwas.nf -c gwas.config -with-docker montilab/gwas
```


### Locally Run Example Data [Centos7]
```bash
$ module load R/3.6.0
$ module load vcftools
$ module load bcftools
$ module load plink/2.00a1LM
$ module load annovar/2018apr

# Modify paths and parameters in gwas.config file

$ nextflow gwas.nf -c gwas.config
```

### Expected Output
```bash

N E X T F L O W  ~  version 19.04.1
Launching `gwas.nf` [jolly_fermi] - revision: 46311ebd05
-

G W A S  ~  P I P E L I N E

================================
indir     : <YOUR PATH>/data/
outdir    : <YOUR PATH>/results

vcf       : <YOUR PATH>/data/toy_vcf.csv
pheno     : <YOUR PATH>/data/pheno_file_logistic.csv
snpset    : <YOUR PATH>/data/snpset.txt

phenotypes: outcome
covars    : age,sex,PC1,PC2,PC3,PC4
model     : logistic
test      : Score
ref       : hg19

-
[warm up] executor > local
executor >  local (142)
[81/6ec45d] process > qc_miss         [100%] 22 of 22 â
[a8/52b952] process > annovar_ref     [100%] 1 of 1 â
[9f/e712d6] process > qc_mono         [100%] 22 of 22 â
[ad/297afe] process > vcf_to_gds      [100%] 22 of 22 â
[fb/0c4c7a] process > merge_gds       [100%] 1 of 1 â
[25/8a55ce] process > pcair           [100%] 1 of 1 â
[d8/8b7b53] process > pcrelate        [100%] 1 of 1 â
[44/5a1f9a] process > caf_by_group    [100%] 22 of 22 â
[b0/28f37d] process > nullmod         [100%] 1 of 1 â
[f1/f8069c] process > gwas            [100%] 22 of 22 â
[de/19afee] process > merge_by_chr    [100%] 22 of 22 â
[08/6939ef] process > combine_results [100%] 1 of 1 â
[14/ac622b] process > annovar_input   [100%] 1 of 1 â
[05/945aff] process > plot            [100%] 1 of 1 â
[2d/76748f] process > annovar         [100%] 1 of 1 â
[28/2ca993] process > add_annovar     [100%] 1 of 1 â
Completed at: 06-Jul-2020 16:07:44
Duration    : 1m 51s
CPU hours   : 0.2
Succeeded   : 142
```
