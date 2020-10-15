#!/usr/bin/env nextflow


/*
USER INPUT PARAMETERS
*/
date = new Date().format( 'yyyyMMdd' )

params.vcf_list    = null
params.pheno       = null

params.phenotype  = null
params.covars      = null

params.pca_grm     = false
params.snpset      = null
params.grm         = null

params.model       = "linear"
params.test        = "Score"

params.gwas        = true
params.imputed     = false

params.gene_based  = false
params.max_maf     = 0.01
params.method      = "Burden"

params.longitudinal = false
params.random_slope = "null"

params.group               = null
params.min_maf             = 0.1
params.max_pval_manhattan  = 0.5

params.max_pval    = 0.01
params.ref_genome  = "hg19"

params.help        = false

println()

/*
OUTPUT DIRECTORY
*/

params.outdir = "$PWD/Analysis_Results-${date}"

/*
HELP MESSAGE
*/

if(params.help){
log.info """
USAGE:

The command lines for runing the pipeline is:

 for gwas:
 nextflow gwas.nf --vcf_list ./data/toy_vcf.csv --pheno ./data/pheno_file_logistic.csv --snpset $PWD/data/snpset.txt --phenotype outcome --covars age,sex,PC1,PC2,PC3,PC4 --ref_genome hg19 --gwas --model linear --test Wald

 for Gene-based Test:

 for longitudinal analysis:

Mandatory arguments:
--vcf_list                 String        Path to the two-column mapping csv file: id , file_path 
--pheno                    String        Path to the phenotype file
--phenotype                String        Name of the phenotype column

Optional arguments:
--outdir                   String        Path to the master folder to store all results
--covars                   String        Name of the covariates to include in analysis model separated by comma (e.g. "age,sex,educ")
--pca_grm                  Logical       If true, run PCAiR (generate PCA in Related individuals) and PCRelate (generate genomic relationship matrix)
--snpset                   String        Path to the two column txt file separated by comma: chr,pos (can only be effective when pca_grm = true)
--grm                      String        Path to the genomic relationship matrix (can only be effective when pca_grm = false)
--model                    String        Name of regression model for gwas: "linear" or "logistic"
--test                     String        Name of statistical test for significance: "Wald" or "Score" ("Wald" test is only effective for linear regression in GWAS)
--gwas                     Logical       If true, run gwas
--imputed                  Logical       If true, use dosages in regression model (DS columns needed in input vcf files)
--gene_based               Logical       If true, run aggregate test for genes based on hg19 reference genome
--max_maf                  Numeric       Threshold for maximun minor allele frequencies of SNPs to be aggregated
--method                   String        Name of aggregation test method: "Burden", "SKAT", "fastSKAT", "SMMAT" or "SKATO"
--longitudinal             Logical       If true, run genome-wide longitudianl analysis
--random_slope             String        if set to "null", random intercept only model is run; else run random slope and random intercept model
--group                    String        Name of the group variable based on which the allele frequencies in each subgroup is calculated (can be left empty)
--min_maf                  Numeric       Threshold for minimun minor allele frequencies of SNPs to include in QQ- and Manhattan-plot
--max_pval_manhattan       Numeric       Threshold for maximun p-value of SNPs to show in Manhattan-plot 
--max_pval                 Numeric       Threshold for maxumun p-value of SNPs to annotate
--ref_genome               String        Name of the reference genome for annotation: hg19 or hg38


"""
  exit 1
}

log.info """\
-

G W A S  ~  P I P E L I N E

================================
outdir    : $params.outdir

vcf       : $params.vcf_list
pheno     : $params.pheno
snpset    : $params.snpset

phenotype : $params.phenotype
covars    : $params.covars
model     : $params.model
test      : $params.test
ref       : $params.ref_genome

-
"""

// Read vcf files
Channel
  .fromPath(params.vcf_list).splitCsv(header: false)
  .map {row -> tuple(row[0], file(row[1]))}
  .ifEmpty {error "File ${params.vcf_list} not parsed properly"}
  .set {vcf_files}

/*
** STEP 0 Configuration Check
*/
if( !(params.pca_grm | params.gwas | params.gene_based | params.longitudinal) ){
  log.info """
  EXIT: NO ANALYSIS SELECTED
  """
  exit 1
}
if(params.gwas&params.gene_based|params.gwas&params.longitudinal|params.gene_based&params.longitudinal|params.gwas&params.gene_based&params.longitudinal){
  log.info """
  EXIT: MORE THAN ONE ANLYSES SELECTED
  """
  exit 1
}

/*
** STEP 1 QC
*/
process qc_miss {
  tag "$chr"
  publishDir "${params.outdir}/QC/qc_miss", mode: 'copy'
  
  input:
  set val(chr), file(vcf) from vcf_files

  output:
  set val(chr), file("*vcf.gz") into qc1

  script:
  """
  vcftools --gzvcf $vcf --max-missing 0.99 --recode --stdout | gzip -c > ${chr}_qc1.vcf.gz
  """
}

process qc_mono {
  tag "$chr"
  publishDir "${params.outdir}/QC/qc_mono", mode: 'copy'
  
  input:
  set val(chr), file(vcf) from qc1

  output:
  set val(chr), file('*vcf.gz') into qc2_1, qc2_2

  script:
  """
  bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $vcf -Oz -o ${chr}_qc2.vcf.gz
  """
}

/*
** STEP 2 vcf_to_gds
*/
process vcf_to_gds {
  tag "$chr"
  publishDir "${params.outdir}/GDS/gds_files", mode: 'copy'
    
  input:
  set val(chr), file(vcf) from qc2_1

  output:
  file '*'
  file('*.gds') into gds_files_1
  set val(chr), file('*.gds') into gds_files_2, gds_files_3

  script:
  """
  02_vcf_to_gds.R $vcf ${chr}.gds ${chr}_vcf_to_gds.log
  """
}

process merge_gds {
  publishDir "${params.outdir}/GDS/gds_merged", mode: 'copy'

  input:
  file(gds_files) from gds_files_1.collect()

  output:
  file '*'
  file 'merged.gds' into gds_merged_1, gds_merged_2, gds_merged_3

  script:
  """
  02_merge_gds.R $gds_files
  """
}

/*
** STEP 3 PCA and GRM
*/
if (params.pca_grm) {
  process pcair {
    publishDir "${params.outdir}/PCA_GRM/pcair", mode: 'copy'
        
    input:
    file(gds_merged) from gds_merged_1

    output:
    file '*' into pcair1, pcair2, pcair3

    script:
    """
    03_PC_AiR.R $gds_merged ${params.pheno} ${params.phenotype} ${params.covars} ${params.snpset}
    """
  }

  process pcrelate {
    publishDir "${params.outdir}/PCA_GRM/pcrelate", mode: 'copy'
    
    input:
    file(gds_merged) from gds_merged_2
    file '*' from pcair1.collect()

    output:
    file '*' into pcrelate

    script:
    """
    03_PC_Relate.R $gds_merged ${params.pheno} ${params.phenotype} ${params.covars} ${params.snpset} analysis.sample.id.rds annot.rds pruned.rds king.rds pcs.rds pc.df.rds
    """
  }
}

/*
** STEP 4 nullmod and gwas/gene-based/longitudinal analysis
*/
if ( (params.gwas | params.gene_based) & params.pca_grm) {
  process nullmod {
    publishDir "${params.outdir}/Association_Test/nullmod", mode: 'copy'
  
    input:
    file(gds_merged) from gds_merged_3
    file '*' from pcair2.collect()
    file '*' from pcrelate.collect()

    output:
    file '*' into nullmod

    script:
    """
    04_nullmod.R $gds_merged ${params.phenotype} ${params.covars} ${params.model} analysis.sample.id.rds annot.rds pc.df.rds grm.rds
    """
  }
}

if ( (params.gwas | params.gene_based) & !params.pca_grm) {
  process nullmod_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/nullmod", mode: 'copy'
  
    input:
    file(gds_merged) from gds_merged_3

    output:
    file '*' into nullmod, nullmod1

    script:
    """
    04_nullmod_skip_pca_grm.R $gds_merged ${params.pheno} ${params.phenotype} ${params.covars} ${params.model} ${params.grm}
    """
  }
}

if (params.gwas & params.pca_grm) {
  process gwas {
    tag "$chr"
    publishDir "${params.outdir}/Association_Test/gwas", mode: 'copy'
    
    input:
    set val(chr), file(gds) from gds_files_2
    file '*' from nullmod.collect()
 
    output:
    set val(chr), file('*.csv') into gwas1
    file '*' into gwas2

    script:
    """
    04_gwas.R ${gds} annot_pc.rds nullmod.rds ${params.test} ${params.imputed} ${chr}.csv ${chr}_gwas.log
    """
  }
}

if (params.gwas & !params.pca_grm) { 
  process gwas_skip_pca_grm {
    tag "$chr"
    publishDir "${params.outdir}/Association_Test/gwas", mode: 'copy'
    
    input:
    set val(chr), file(gds) from gds_files_2
    file '*' from nullmod.collect()
 
    output:
    set val(chr), file('*.csv') into gwas1
    file '*' into gwas2

    script:
    """
    04_gwas.R ${gds} annot.rds nullmod.rds ${params.test} ${params.imputed} ${chr}.csv ${chr}_gwas.log
    """
  }
}

if (params.gene_based & params.pca_grm) {
  process gene_based {
    tag "$chr"
    publishDir "${params.outdir}/Association_Test/gene_based", mode: 'copy'
    
    input:
    set val(chr), file(gds) from gds_files_2
    file '*' from nullmod.collect()
 
    output:
    set val(chr), file('*.csv') into gene_based
    file '*'
  
    script:
    """
    04_gene_based.R ${gds} annot_pc.rds nullmod.rds ${params.max_maf} ${params.method} ${chr}.csv ${chr}.rds ${chr}_gene_based.log
    """
  }
}

if (params.gene_based & !params.pca_grm) {
  process gene_based_skip_pca_grm {
    tag "$chr"
    publishDir "${params.outdir}/Association_Test/gene_based", mode: 'copy'
    
    input:
    set val(chr), file(gds) from gds_files_2
    file '*' from nullmod.collect()
 
    output:
    set val(chr), file('*.csv') into gene_based
    file '*'

    script:
    """
    04_gene_based.R ${gds} annot.rds nullmod.rds ${params.max_maf} ${params.method} ${chr}.csv ${chr}.rds ${chr}_gene_based.log
    """
  }
}

if (params.longitudinal&params.pca_grm) {
  process nullmod_longitudinal {
    publishDir "${params.outdir}/Association_Test/nullmod_longitudinal", mode: 'copy'
  
    input:
    file '*' from pcair2.collect()
    file '*' from pcrelate.collect()

    output:
    file '*' into nullmod_longitudinal

    script:
    """
    04_nullmod_longitudinal.R ${params.pheno} ${params.phenotype} ${params.covars} ${params.model} analysis.sample.id.rds pc.df.rds grm.rds ${params.random_slope}
    """
  }
}

if (params.longitudinal & !params.pca_grm) {
  process nullmod_longitudinal_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/nullmod_longitudinal", mode: 'copy'
  
    input:
    file(gds_merged) from gds_merged_3

    output:
    file '*' into nullmod_longitudinal
    file '*' into nullmod1

    script:
    """
    04_nullmod_longitudinal_skip_pca_grm.R $gds_merged ${params.pheno} ${params.phenotype} ${params.covars} ${params.model} ${params.grm} ${params.random_slope}
    """
  }
}

if (params.longitudinal) {
  process gwas_longitudinal {
    tag "$chr"
    publishDir "${params.outdir}/Association_Test/gwla", mode: 'copy'
    
    input:
    set val(chr), file(gds) from gds_files_2
    file '*' from nullmod_longitudinal.collect()
 
    output:
    set val(chr), file('*.txt') into gwas_longitudinal
   
    script:
    """
    04_gwas_longitudinal.R ${gds} nullmod_longitudinal.rds ${chr}.txt ${chr}_gwas_longitudinal.log
    """
  }
}

/*
** STEP 5 summary and plot
*/

if (params.gene_based) {
  process combine_results_gene {
    publishDir "${params.outdir}/Summary_Plot/combined_results", mode: 'copy'
  
    input:
    file(gene) from gene_based.collect()
  
    output:
    file '*' into combined_results

    script:
    """
    05_combine_results_gene.R ${gene}
    """
  }

  process plot_gene {
    publishDir "${params.outdir}/Summary_Plot/qq_plot", mode: 'copy'
  
    input:
    file '*' from combined_results.collect()
  
    output:
    file '*' into qq_plot

    script:
    """
    05_qqplot_gene.R all_chr.csv
    """
  }
}

if ( (params.gwas | params.longitudinal) & params.pca_grm) {
 process caf_by_group {
  tag "$chr"
  publishDir "${params.outdir}/Summary_Plot/caf_by_group", mode: 'copy'
  
  input:
  set val(chr), file(gds) from gds_files_3
  file '*' from pcair3.collect()
  
  output:
  file '*' 
  set val(chr), file('*.csv') into caf_by_group

  script:
  """
  05_caf_by_group.R ${gds} ${params.pheno} analysis.sample.id.rds ${params.model} ${params.phenotype} ${params.group} ${chr}_caf_by_group.csv ${chr}_caf_by_group.log
  """
 }
}

if ( (params.gwas | params.longitudinal) & !params.pca_grm) {
 process caf_by_group_skip_pca_grm {
  tag "$chr"
  publishDir "${params.outdir}/Summary_Plot/caf_by_group", mode: 'copy'
  
  input:
  set val(chr), file(gds) from gds_files_3
  file '*' from nullmod1.collect()
  
  output:
  file '*' 
  set val(chr), file('*.csv') into caf_by_group

  script:
  """
  05_caf_by_group.R ${gds} ${params.pheno} analysis.sample.id.rds ${params.model} ${params.phenotype} ${params.group} ${chr}_caf_by_group.csv ${chr}_caf_by_group.log
  """
 }
}

if (params.longitudinal) {
  process coincide_gwas {
    tag "$chr"
    publishDir "${params.outdir}/Summary_Plot/gwla", mode: 'copy'
  
    input:
    set val(chr), file(txt) from gwas_longitudinal
  
    output:
    file '*' 
    set val(chr), file('*.csv') into gwas1

    script:
    """
    05_coincide_by_chr.R ${txt} ${chr}.csv ${chr}_coincide.log
    """
  }
}

if (params.gwas | params.longitudinal) {
  process merge_by_chr {
    tag "$chr"
    publishDir "${params.outdir}/Summary_Plot/merge_by_chr", mode: 'copy'
  
    input:
    set val(chr), file(csv1) from gwas1
    set val(chr), file(csv2) from caf_by_group
  
    output:
    file '*' 
    file('*.csv') into merge_by_chr

    script:
    """
    05_merge_by_chr.R ${csv1} ${csv2} ${params.model} ${chr}_caf_annotated.csv ${chr}_merge.log
    """
  }
}

if(params.gwas|params.longitudinal){
process combine_results {
  publishDir "${params.outdir}/Summary_Plot/combined_results", mode: 'copy'
  
  input:
  file(csv_files) from merge_by_chr.collect()
  
  output:
  file '*' into combined_results1
  file '*' into combined_results2

  script:
  """
  05_combine_results.R ${csv_files}
  """
}

process plot {
  publishDir "${params.outdir}/Summary_Plot/qq_manhattan", mode: 'copy'
  
  input:
  file '*' from combined_results1.collect()
  
  output:
  file '*' into qq_manhattan

  script:
  """
  05_qqplot_manhattanplot.R all_chr_caf_annotated.csv ${params.min_maf} ${params.max_pval_manhattan}
  """
}

/*
** STEP 6 annotation
*/

process annovar_input {
  publishDir "${params.outdir}/Annotation/annovar_input", mode: 'copy'
    
  input:
  file '*' from combined_results2.collect()
  
  output:
  file '*' into annovar_input1
  file '*' into annovar_input2

  script:
  """
  06_annovar_input.R ${params.max_pval}
  """
}

process annovar_ref {  
  publishDir "${params.outdir}/Annotation/", mode: 'copy'
  
  script:
  """
  mkdir -p ${params.outdir}/Annotation
  annotate_variation.pl --downdb --buildver ${params.ref_genome} --webfrom annovar refGene ${params.outdir}/Annotation/humandb
  """
}

process annovar {
  publishDir "${params.outdir}/Annotation/annovar", mode: 'copy'
    
  input:
  file '*' from annovar_input1.collect()

  output:
  file '*' into annovar

  script:
  """
  table_annovar.pl -build ${params.ref_genome} top_snps_input.txt ${params.outdir}/Annotation/humandb/ -out top_annotation -remove -protocol refGene -operation g -nastring . -csvout
  """
}

process add_annovar {
  publishDir "${params.outdir}/Annotation/annotated_results", mode: 'copy'
  
  input:
  file '*' from annovar_input2.collect()
  file '*' from annovar.collect()

  output:
  file '*' into report

  script:
  """
  06_add_anno_results.R ${params.ref_genome}
  """
}
}

/*
** STEP 7 Report
*/

if(params.gwas|params.longitudinal){
  process report {
    publishDir "${params.outdir}/Report", mode: 'copy'
  
    input:
    file '*' from report.collect()
  
    output:
    
    script:
    """
    mkdir -p ${params.outdir}/Report
    Rscript -e 'rmarkdown::render("$PWD/bin/07_report.Rmd")' ${params.outdir}
    mv $PWD/bin/07_report.html ${params.outdir}/Report/Report.html
    """
  }
}

if(params.gene_based){
  process report_gene {
    publishDir "${params.outdir}/Report", mode: 'copy'
  
    input:
    file '*' from qq_plot.collect()
  
    output:

    script:
    """
    mkdir -p ${params.outdir}/Report
    Rscript -e 'rmarkdown::render("$PWD/bin/07_report_gene.Rmd")' ${params.outdir} ${params.max_pval}
    mv $PWD/bin/07_report_gene.html ${params.outdir}/Report/Report.html
    """
  }
}

