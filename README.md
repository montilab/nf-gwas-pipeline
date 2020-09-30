# nf-gwas-pipeline
A Nextflow Genome-Wide Association Study (GWAS) Pipeline

### Clone Repository
```bash

$ git clone https://github.com/gurinovich/GWASpipeline

```
### Initalize Paths to Test Data
```bash

$ cd GWASpipeline/pipeline 
$ python utils/paths.py  

```

### Download Nextflow Executable
```

$ curl -s https://get.nextflow.io | bash

```

### Locally Run Example Data [Centos7]
```bash
$ module load R/3.6.0
$ module load vcftools
$ module load bcftools
$ module load plink/2.00a1LM
$ module load annovar/2018apr
$ mkdir results

# Modify paths and parameters in gwas.config file

$ nextflow gwas.nf -c gwas.config
```

### Run with docker
```
$ docker pull montilab/gwas:latest

$ nextflow gwas.nf -c gwas.config -with-docker montilab/gwas
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

vcf       : <YOUR PATH>/data//toy_vcf.csv
pheno     : <YOUR PATH>/data//pheno_file_logistic.csv
snpset    : <YOUR PATH>/data//snpset.txt

phenotypes: outcome
covars    : age,sex,PC1,PC2,PC3,PC4
model     : logistic
test      : Score
ref       : hg19

-
[warm up] executor > local
executor >  local (142)
[81/6ec45d] process > qc_miss         [100%] 22 of 22 ✔
[a8/52b952] process > annovar_ref     [100%] 1 of 1 ✔
[9f/e712d6] process > qc_mono         [100%] 22 of 22 ✔
[ad/297afe] process > vcf_to_gds      [100%] 22 of 22 ✔
[fb/0c4c7a] process > merge_gds       [100%] 1 of 1 ✔
[25/8a55ce] process > pcair           [100%] 1 of 1 ✔
[d8/8b7b53] process > pcrelate        [100%] 1 of 1 ✔
[44/5a1f9a] process > caf_by_group    [100%] 22 of 22 ✔
[b0/28f37d] process > nullmod         [100%] 1 of 1 ✔
[f1/f8069c] process > gwas            [100%] 22 of 22 ✔
[de/19afee] process > merge_by_chr    [100%] 22 of 22 ✔
[08/6939ef] process > combine_results [100%] 1 of 1 ✔
[14/ac622b] process > annovar_input   [100%] 1 of 1 ✔
[05/945aff] process > plot            [100%] 1 of 1 ✔
[2d/76748f] process > annovar         [100%] 1 of 1 ✔
[28/2ca993] process > add_annovar     [100%] 1 of 1 ✔
Completed at: 06-Jul-2020 16:07:44
Duration    : 1m 51s
CPU hours   : 0.2
Succeeded   : 142

```
