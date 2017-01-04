from __future__ import division
from operator import add
from math import log
import csv
import pickle
import os
import sys
import glob


def makeGenotype(line,idCol):
    AA=map(float,line[5::3])
    AB=map(float, line[6::3])
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return (line[idCol], genotype)


def makeGenotypeCheckRef(line, checkMap):
    rsid=line[0]
    gen=line[1]
    try:
        if checkMap[rsid]=="keep":
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))


        elif checkMap[rsid]=="flip":
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))
        else:
            genotype=[0.0]*(len(gen)/3)

    except KeyError:
        print("SNP {} was not accounted for in the alignment checking step, discarding this SNP".format(rsid))
        genotype=[0.0]*(len(gen)/3)
    finally:
        return (rsid, genotype)
def getMaf(geno):
    AA=geno[0::3]
    AB=geno[1::3]
    AA2=[x*2 for x in AA]
    A1count=map(add, AA2, AB)
    A1F=sum(A1count)/(float(len(AA2))*2)
    return A1F


def filterGWASByP(GWASRdd, pcolumn,  pHigh, oddscolumn,idcolumn,  logOdds,pLow=0):
    GWAS_Pfiltered=GWASRdd.filter(lambda line: (float(eval(line[pcolumn]))<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair

def filterGWASByP_DF(GWASdf, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASdf.rdd.filter(lambda line: (float(line[pcolumn])<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        print("Taking the log of odds ratios")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:
        print("Using the original values of effect sizes")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair


def multiplyOdds(genotypeRDDLine, oddsMap2):
    if genotypeRDDLine is not None:
        if genotypeRDDLine[0] in oddsMap2:
            OR=oddsMap2[genotypeRDDLine[0]]
            return [x*OR for x in genotypeRDDLine[1]]
        else:
            pass
    else:
        pass


def checkAlignmentDF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1]
    genoA2=dataframe[2]
    genoA1F=dataframe[3]
    gwasA1=dataframe[5]
    gwasA2=dataframe[6]
    gwasA1F=dataframe[7]
    try:
        if genoA1==bpMap[genoA1]:
            if gwasA1F==".":
                flag="discard"
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5

                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff*genoA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"
        print("Invalid Genotypes for SNP {}".format(snpid))

    finally:
        return (snpid,flag)



def makeGenotype(line):
    AA=map(float,line[5::3])
    AB=map(float, line[6::3])
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return (line[1], genotype)

def makeGenotypeCheckRef(line, checkMap, toDF=False):
    rsid=line[0]
    gen=line[1]
    try:
        if checkMap[rsid]=="keep":
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))


        elif checkMap[rsid]=="flip":
            AA=gen[0::3]
            AB=gen[1::3]
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))
        else:
            genotype=[0.0]*(len(gen)/3)

    except KeyError:
        print("SNP {} was not accounted for in the alignment checking step, discarding this SNP".format(rsid))
        genotype=[0.0]*(len(gen)/3)
    finally:
        if toDF:
            return [rsid]+genotype
        else:
            return (rsid, genotype)


def filterGWASByP(GWASRdd, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASRdd.filter(lambda line: (float(eval(line[pcolumn]))<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        print("Filtering GWAS, taking the log of odds ratios")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:
        print("Filtering GWAS, keeping the original values of effect sizes")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair

def filterGWASByP_DF(GWASdf, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASdf.rdd.filter(lambda line: (float(line[pcolumn])<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        print("Filtering GWAS, taking the log of odds ratios")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],log(float(line[oddscolumn]))))
    else:
        print("Filtering GWAS, keeping the original values of effect sizes")
        GWAS_Odds=GWAS_Pfiltered.map(lambda line: (line[idcolumn],  float(line[oddscolumn])))
    GWASoddspair=GWAS_Odds.collectAsMap()
    return GWASoddspair


def multiplyOdds(genotypeRDDLine, oddsMap2):
    if genotypeRDDLine is not None:
        if genotypeRDDLine[0] in oddsMap2:
            OR=oddsMap2[genotypeRDDLine[0]]
            return [x*OR for x in genotypeRDDLine[1]]
        else:
            pass
    else:
        pass

def checkAlignmentDF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1]
    genoA2=dataframe[2]
    genoA1F=dataframe[3]
    gwasA1=dataframe[5]
    gwasA2=dataframe[6]
    gwasA1F=dataframe[7]
    try:
        if genoA1==bpMap[genoA1]:
            if gwasA1F==".":
                flag="discard"
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5

                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff*genoA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"
        print("Invalid Genotypes for SNP {}".format(snpid))

    finally:
        return (snpid,flag)

def checkAlignmentDFnoMAF(dataframe, bpMap):
    snpid=dataframe[0]
    genoA1=dataframe[1]
    genoA2=dataframe[2]
    gwasA1=dataframe[4]
    gwasA2=dataframe[5]

    try:
        if genoA1==bpMap[genoA1]:
            flag="discard"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"
        print("Invalid Genotypes for SNP {}".format(snpid))

    finally:
        return (snpid,flag)





def checkAlignment(line):
    bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
    genoA1=line[0][0][0]
    genoA2=line[0][0][1]
    genoA1F=line[0][1]
    gwasA1=line[1][0][0]
    gwasA2=line[1][0][1]
    gwasA1F=line[1][1]

    try:
        if genoA1==bpMap[genoA1]:
            if gwasA1F==".":
                flag="discard"
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5

                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff*genoA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"

    finally:
        return flag



def getSampleNames(sampleFileName, sampleDelim, sampleIDCol, skip=0):
    labels=[]
    with open(sampleFileName, "r") as f:
        subjList=[item.split(sampleDelim) for item in f.read().splitlines()]
        for i in range(len(sampleIDCol)):
            subjNames=[x[i] for x in subjList[skip::]]
            subjNames=[name.strip('"') for name in subjNames]
            column=["Label"+str(i)]+subjNames
            labels.append(column)
    return labels




# output each PRS for each sample, in the form of [sample, *scores],
# and a list of pvalues that are in the order of the scores from each p-value


def writePRS(prsResults, outputFile, samplenames=None, dialect=None):
    samplesize=len(prsResults.values()[0][1])
    if not samplenames:
        print "No sample names provided, generating sample names"
        samplenames={"Label":["Sample"+str(i) for i in range(len(prsResults))]}

    outputdata=samplenames
    for pvalue in prsResults.keys():
        outputdata.append(["SNP_count_{}".format(pvalue)]+[prsResults[pvalue][0]]*samplesize)
        outputdata.append(["PRS_{}".format(pvalue)]+prsResults[pvalue][1])

    try:
        with open(outputFile, "w") as f:
            csvwriter=csv.writer(f, dialect=dialect)
            for row in zip(*outputdata):
                csvwriter.writerow(row)
            print("Successfully wrote scores to "+ os.path.basename(outputFile))
    except:
        e = sys.exc_info()[0]
        print( "<p>Error: %s</p>" % e )
        print("Data output was unsuccessful.")
        print("All is not lost, final results saved as binary format in file 'PRSOutput.pk'")
        with open(os.path.dirname(outputFile)+"/PRSOutput.pk", "wb") as f:
            pickle.dump(prsResults, f)
    return outputdata

## The following function is for checking allignment in genotype RDDs
def checkAlignment(line):
    bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
    genoA1=line[0][0][0]
    genoA2=line[0][0][1]
    genoA1F=line[0][1]
    gwasA1=line[1][0][0]
    gwasA2=line[1][0][1]
    gwasA1F=line[1][1]

    try:
        if genoA1==bpMap[genoA1]:
            if gwasA1F==".":
                flag="discard"
            else:
                gwasA1F=float(gwasA1F)
                genoA1Fdiff=genoA1F-0.5
                gwasA1Fdiff=float(gwasA1F)-0.5

                if abs(genoA1Fdiff)<0.1 or abs(gwasA1Fdiff)<0.1:
                    flag="discard"
                else:
                    if genoA1Fdiff*genoA1Fdiff>0:
                        flag="keep"
                    else:
                        flag="flip"
        elif genoA2==gwasA1 or genoA2==bpMap[gwasA1]:
            flag="flip"
        else:
            flag="keep"

    except KeyError:
        flag="discard"

    finally:
        return flag





# output each PRS for each sample, in the form of [sample, *scores],
# and a list of pvalues that are in the order of the scores from each p-value
