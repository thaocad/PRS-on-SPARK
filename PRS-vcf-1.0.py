# coding: utf-8

# ### Calculating PRS using VCF files

from __future__ import division
from pyspark import SparkConf, SparkContext
from operator import add
import re
import glob, os
import csv
from collections import Counter
import ntpath
import functools
#from functools import reduce
from math import log
import itertools
import PRS_VCF_utils
APP_NAME="MAVANvcfPRS"


# PRS calculation on pruned NFP data, 50, 5, VIF=2

#**ATTN: python index starts at 0, so if you want to specify the second column, use 1
#**ATTN: please remove the header of the GWAS file if there is any

# define column number for contents in GWAS

gwas_id=0    # column of SNP ID
gwas_p=7     # column of P value
gwas_or=5    # column of odds ratio
gwas_a1=3    # column of a1 in the GWAS
gwas_maf= 10 # column index of maf in the GWAS

# defin column number for contents in genfile
geno_id=2  # column number with rsID
geno_start=9 # column number of the 1st genotype, in the raw vcf files, after separated by the delimiter of choice
geno_a1 = 3  # column number that contains the reference allele

# List of thresholds:
thresholds=[0.5, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001]

# file delimiters:
GWAS_delim="\t"
GENO_delim="\t"

# file names:
home="/home/nyao111/MAVAN_imputed_161121/MOMS_info03_vcf/"  #define homefolder path

gwasFiles="/home/nyao111/PRS_imputed/pgc.mdd.clump.withAF.txt"       # Name of GWAS file


def getFileFromPattern(*pattern): # Multiple patterns need to be put into list format
    files=[]
    for pathpattern in pattern:
        files=glob.glob(files)

genoFileNamePattern=home+"*info03.vcf"

genoFileNames=glob.glob(genoFileNamePattern)
# Alternatively, directly specify filename:
#genoFileName=[home+"fcgene_out_chr21comb.bierut1M_plus_filtered_chr21_c1_EA_COGA.gen",
              #home+"fcgene_out_chr21comb.bierut1M_plus_filtered_chr21_c1_EA_COGEND.gen",
              #home+"fcgene_out_chr22comb.bierut1M_plus_filtered_chr22_c1_EA_COGA.gen",
              #home+"fcgene_out_chr22comb.bierut1M_plus_filtered_chr22_c1_EA_COGEND.gen"]

genoExtension=".vcf"


# programme parameters
log_or=True  # sepcify whether you want to log your odds ratios
check_ref=True # if you know that there are mismatch between the top strand in the genotypes and that of the GWAS, set True. Not checking the reference allele will improve the speed

# sample file path and name
sampleFilePath=home+"MAVAN_35_impute161121_MOM_orderedSamples.csv" # include the full/relative path and name of the sample file
sampleFileDelim=","  # sample File Delimiter
sampleFileID=0   # which column in the sample file has the ID
sample_skip=1  # how many lines to skip so that the sample names can be matched to the genotypes 1-to-1, taking into account the header of the sample file
##output file information

outputPath=home+"MAVAN_MOMS_mddPRS_161125.csv"



# In[2]:

import pyspark
from pyspark.sql import SQLContext

# We can give a name to our app (to find it in Spark WebUI) and configure execution mode

app_name = "PRS"

conf = pyspark.SparkConf().setAppName(app_name)
sc = pyspark.SparkContext(conf=conf)
print(sc)
sc.setLogLevel("WARN")
log4jLogger = sc._jvm.org.apache.log4j
LOGGER = log4jLogger.LogManager.getLogger(__name__)
LOGGER.info("Start Reading Files")



print("="*40)
print("Using these genoytpe files: ")

counter = 0
for filename in genoFileNames:
    if counter<20:
        counter+=1
        print(filename)
    else:
        print("and more....")
        break



genodata=sc.textFile(genoFileNamePattern)
gwasfile=sc.textFile(gwasFiles)
print("Using the GWAS file: {}".format(ntpath.basename(gwasFiles)))
gwastable=gwasfile.filter(lambda line: "snpid" not in line).map(lambda line: line.split(GWAS_delim))
gwastableCA=gwastable.cache()


# ### 1.1 Filter GWAS and prepare odds ratio
#

# In[5]:

maxThreshold=max(thresholds)
gwasOddsMapMax=PRS_VCF_utils.filterGWASByP(GWASRdd=gwastableCA, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=maxThreshold, logOdds=log_or)
gwasOddsMapMaxCA=sc.broadcast(gwasOddsMapMax).value


# ### 2. Initial processing

# In[6]:

# at this step, the genotpe is already filtered to keep only the ones in 'gwasOddsMap'
genointermediate=genodata.filter(lambda line: ("#" not in line)).map(lambda line: line.split(GENO_delim)).filter(lambda line: line[geno_id] in gwasOddsMapMaxCA).map(lambda line: line[0:5]+[chunk.split(":")[3] for chunk in line[geno_start::]]).map(lambda line: line[0:5]+[triplet.split(",") for triplet in line[5::]])

genoAlleles=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1])))
genotable=genointermediate.map(lambda line: (line[geno_id], list(itertools.chain.from_iterable(line[5::])))).mapValues(lambda geno: [float(x) for x in geno])


# In[7]:

print genointermediate.count()


# ### 2.1 Calculate and store MAF

# In[8]:

reload(PRS_VCF_utils)
genoa1f=genointermediate.map(lambda line: (line[geno_id], (line[geno_a1], line[geno_a1+1]), [float(x) for x in list(itertools.chain.from_iterable(line[5::]))])).map(lambda line: (line[0], line[1],PRS_VCF_utils.getMaf(line[2])))

#genoa1f.map(lambda line:"\t".join([line[0], "\t".join(line[1]), str(line[2])])).saveAsTextFile("../MOMS_info03_maf")
genoa1f.first()


# ### 3. Determine whether each SNP needs to be flipped

# In[9]:

checktable=genoa1f.map(lambda line: (line[0], (line[1], line[2]))).join(gwastableCA.map(lambda line:(line[gwas_id], ((line[gwas_a1], line[gwas_a1+1]), line[gwas_maf]))))


# In[10]:

get_ipython().magic(u'time')
reload(PRS_VCF_utils)
flagMap=checktable.mapValues(lambda line: PRS_VCF_utils.checkAlignment(line)).collectAsMap()


# In[15]:

# filter the raw genotype file
reload(PRS_VCF_utils)
if check_ref:
    print("Calculating genotype dosage while taking into account strand alignment differences")
    genotypeMax=genotable.map(lambda line: PRS_VCF_utils.makeGenotypeCheckRef(line, checkMap=flagMap)).cache()
    samplesize=len(genotypeMax.first()[1])
else:
    genotypeMax=genotable.map(lambda line: PRS_VCF_utils.makeGenotype(line, gwasOddsMapCA)).cache()
    samplesize=len(genotypeMax.first()[1])


# In[16]:

# Calculate the PRS with the maximum threshold
# calculate PRS from genotype
reload(PRS_VCF_utils)
def calcPRSFromGeno(genotypeRDD, oddsMap):
    totalcount=genotypeRDD.count()
    multiplied=genotypeRDD.map(lambda line:[call * oddsMap[line[0]] for call in line[1]])
    filtered=multiplied.filter(lambda line: line is not None)
    PRS=multiplied.reduce(lambda a,b: map(add, a, b))
    normalizedPRS=[x/totalcount for x in PRS]
    return (totalcount,PRS)

prsMax=calcPRSFromGeno(genotypeMax, gwasOddsMapMaxCA)
prsDict={}
prsDict[maxThreshold]=prsMax
# Calculate PRS for the rest of the thresholds


# In[29]:

from time import time

def calcNoMax(genotypeRDD, gwasRDD, thresholdlist, prsMap):

    if len(thresholdlist)>1:
        thresholdListNoMax=[x for x in thresholds if x != maxThreshold]
        thresholdNoMaxSorted=sorted(thresholdListNoMax, reverse=True)
    else:
        thresholdNoMaxSorted=thresholdlist
    start=time()
    for threshold in thresholdNoMaxSorted:
        tic=time()
        gwasFiltered=PRS_VCF_utils.filterGWASByP(GWASRdd=gwasRDD, pcolumn=gwas_p, idcolumn=gwas_id, oddscolumn=gwas_or, pHigh=threshold, logOdds=log_or)
        print("filtering GWAS at threshold of {} took {:3.2f} seconds".format( str(threshold), time()-tic) )

        checkpoint=time()-start

        gwasFilteredBC=sc.broadcast(gwasFiltered)
        filteredgenotype=genotypeRDD.filter(lambda line: line[0] in gwasFilteredBC.value)
        if not filteredgenotype.isEmpty():
            prsOther=calcPRSFromGeno(filteredgenotype, gwasFilteredBC.value)
            prsMap[threshold]=prsOther
            print("finished calculating PRS at threshold of {}, used {:3.2f} seconds".format(str(threshold), time()-checkpoint))

    return prsMap

finalresult=calcNoMax(genotypeMax,gwastableCA, thresholds, prsDict)


# In[34]:

reload(PRS_VCF_utils)
subjNames=PRS_VCF_utils.getSampleNames(sampleFilePath,sampleFileDelim,sampleFileID, skip=1)


# In[60]:

os.path.dirname(outputPath)


# In[66]:

reload(PRS_VCF_utils)
test=PRS_VCF_utils.writePRS(finalresult,  outputPath, samplenames=subjNames)
