# PRS-on-SPARK
Polygenic risk score pipeline for large genotype data, generally imputed data

## Installation

To clone the repository, use 
```
git clone https://github.com/seriousNel/PRS-on-SPARK.git
```

## Software Requirements

The notebooks and scripts require the following to run :
+ spark-2.0.0 +
+ Pyhon 2.7

Instruction of installing Apache Spark can be found [here](https://www.santoshsrinivas.com/installing-apache-spark-on-ubuntu-16-04/)


## What the pipeline does:
+ Calculate PRS from a genotype file (in .gen or .vcf format) and a GWAS file
+ Correct the strand alignment descrpencies between genotype  and GWAS data. 

## What it cannot do :
+ Performs quality control of genotype data

## Default format requirement:
### GWAS
By default, the GWAS should have the same format as that of a GWAS file obtained from Psychiatric Genomics Consortium (PGC). 

|     snpid|hg18chr|      bp| a1| a2|    or|    se|  pval| info|ngt|    CEUaf |
|----------|-------|--------|---|---|------|------|------|-----|---|----------|
| rs2241028|      2|74612442|  A|  G|0.9991|0.0255|0.9723|0.968|  0| 0.0688073|
| rs6707302|      2|74614930|  T|  C|0.9762|0.0179|0.1779|0.931|  0|   0.87963|
|rs13428246|      2|74619088|  T|  C|0.9945|0.2186|0.9799|0.405|  0|0.00458716|
|  rs715407|      2|74619201|  T|  G| 1.025| 0.018|0.1734|0.947|  0|  0.119266|
|rs12104595|      2|74621725|  T|  C|0.9777|0.1524|0.8823|0.918|  0|0.00458716|

You can change your GWAS to the same format, or use optional parameter flags to let the script know about the format you are using. More details below.

### .gen file
from [www.shapeit.fr](http://www.shapeit.fr/pages/m02_formats/gensample.html) :
A .gen file is a SPACE delimited file. Each line corresponds to a single SNP. The first 5 columns are:
Chromosome number [integer]
SNP ID [string]
SNP physical position (bp) [integer]
First allele [string]
Second allele [string]

### .vcf file 
This is a default format for the genotype data returned from Sanger Institute. Details about the format can be found [here](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) 

## Running command-line script PRS_run.py
### Parameters

There are three required positional parameters: 
1. path to the genotype file(s)
2. path to the GWAS file
3. path to the output file


A description of the parameters for the script can be obtained by typing: 
```
python PRS_run.py
```
which gives: 

```

usage: PRS_run.py [-h] [--gwas_id GWAS_ID] [--gwas_p GWAS_P]
                  [--gwas_or GWAS_OR] [--gwas_a1 GWAS_A1]
                  [--gwas_maf GWAS_MAF] [--filetype {GEN,VCF}]
                  [--thresholds THRESHOLDS [THRESHOLDS ...]]
                  [--GWAS_delim GWAS_DELIM] [--GWAS_no_header] [--log_or]
                  [--check_ref] [--app_name APP_NAME]
                  [--sample_file SAMPLE_FILE] [--sample_delim SAMPLE_DELIM]
                  [--sample_file_ID SAMPLE_FILE_ID [SAMPLE_FILE_ID ...]]
                  [--sample_file_skip SAMPLE_SKIP] [--no_maf]
                  GENO GWAS Output

PRS Script Parameters

positional arguments:
  GENO                  Name of the Genotype files, can be a name or path, or
                        name patterns with '*'
  GWAS                  Name of the GWAS file, can be a name or path.
  Output                The path and name for the output file

optional arguments:
  -h, --help            show this help message and exit
  --gwas_id GWAS_ID     Column number in your GWAS that contains SNP ID, with
                        first column being 0, default is 0
  --gwas_p GWAS_P       Column number in your GWAS that contains p-value, with
                        first column being 0, default is 7
  --gwas_or GWAS_OR     Column number in your GWAS that contains odds-
                        ratio/beta, with first column being 0, default is 5
  --gwas_a1 GWAS_A1     Column number in your GWAS that contains allele A1,
                        with first column being 0, default is 3. Allele A2 is
                        assumed to be at column [gwas_a1+1]
  --gwas_maf GWAS_MAF   Column number in your GWAS that contains frequency of
                        A1, with first column being 0, default is 10.
  --filetype {GEN,VCF}  The type of genotype file used as inputm choose
                        between VCF and GEN, default is VCF
  --thresholds THRESHOLDS [THRESHOLDS ...]
                        The p-value thresholds that controls which SNPs are
                        used from the GWAS. Specifying the p-values simply by
                        input one after another. default is [0.5, 0.2, 0.1,
                        0.05, 0.01, 0.001, 0.0001]
  --GWAS_delim GWAS_DELIM
                        Delimtier of the GWAS file, default is tab-delimiter
  --GWAS_no_header      Adding this parameter signals that there is no headers
                        for the GWAS. The default is to assume that GWAS has
                        column names
  --log_or              Adding this parameter tells the script to log the
                        effect sizes provided in the GWAS
  --check_ref           Adding this option tells the script to theck reference
                        allele when determining genoypte calls. Default is not
                        checking
  --app_name APP_NAME   Give your spark application a name. Default is PRS.
  --sample_file SAMPLE_FILE
                        path and name of the file that contain the sample
                        labels. It is assumed that the sample labels are
                        already in the same order as in the genotype file.
  --sample_delim SAMPLE_DELIM
                        Delimiter of the sample file. Default is comma
  --sample_file_ID SAMPLE_FILE_ID [SAMPLE_FILE_ID ...]
                        Specify which columns in the sample file are used as
                        labels. Can use one integer to specify one column, or
                        multiple integers to specify multiple columns. Default
                        is the first column
  --sample_file_skip SAMPLE_SKIP
                        Specify how many lines to skip in the sample file,
                        i.e. which row do the labels start. Default is 1,
                        which assumes that the sample files has column names
                        and the labels start on the second line
  --no_maf              By default, the pipeline calculated the allele
                        frequency in the genotype population. Use this flag to
                        tell the script NOT to calculate MAF in the provided
                        propulation and compare it with MAF in the GWAS, e.g,
                        when the GWAS does not provide information for allele
                        frequencies. MAF is needed to check the reference
                        alleles of ambiguous SNPs (those whose A1 and A2 are
                        reverese complements). Not using this will result in
                        ambiguous SNPs be discarded.


```

### Examples:
To calculate PRS from a series of .vcf files, while checking the allele allignment between the genotype and the GWAS, and take the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:
```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --sample_file samplefile.csv --sample_file_id 0 --check_ref --log_or --thresholds  0.2 0.1 0.05
```
To calculate PRS from a series of .gen files, without checking allele alignments, using a GWAS with no header, and taking the log of effect sizes, using p-value thresholds of 0.2, 0.1 , 0.05:

```
spark-submit PRS_run.py "VCF_number*.vcf" pgc.mdd.clump.txt output.csv --filetype GEN --sample_file samplefile.csv --sample_file_id 0 --GWAS_no_header --log_or --thresholds  0.2 0.1 0.05
```

