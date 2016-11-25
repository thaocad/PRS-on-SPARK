from __future__ import division
from operator import add
from math import log

def makeGenotype(line,idCol):
    AA=map(float,line[5::3])
    AB=map(float, line[6::3])
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return (line[idCol], genotype)

bpPair={"C":"G", "G":"C", "A":"T", "T":"A"}

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



def filterGWASByP(GWASRdd, pcolumn,  pHigh, oddscolumn,idcolumn, pLow=0, logOdds=False):
    GWAS_Pfiltered=GWASRdd.filter(lambda line: (float(line[pcolumn])<=pHigh) and (float(line[pcolumn])>=pLow))
    if logOdds:
        print("Taking the log of odds-ratios")
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

def checkAlignment(line):
    bpMap={"A":"T", "T":"A", "C":"G", "G":"C"}
    genoA1=line[0][0][0]
    genoA2=line[0][0][1]
    genoA1F=line[0][1]
    gwasA1=line[1][0][0]
    gwasA2=line[1][0][1]
    gwasA1F=line[1][1]
    flag="keep"
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
        elif genoA1==gwasA2 or genoA1==bpMap[gwasA2]:
            flag="flip"

    except KeyError:
        flag="discard"

    return flag



def getSampleNames(scores, sampleFileName, sampleDelim, sampleIDCol, skip=True):
    samplesize=len(scores)
    with open(sampleFileName, "r") as f:
        subjList=[item.split(sampleDelim) for item in f.read().splitlines()]
        subjNames=[x[sampleIDCol] for x in subjList[skip::]]

    if len(subjNames) == samplesize:
        return subjNames
    else:
        print("Number of Subjects does not equal to number of scores. Using self-made labels")
        selfLabels = ["Subj"+str(i+1) for i in range(samplesize)]
        return selfLabels

def getMaf(geno):
    AA=geno[0::3]
    AB=geno[1::3]
    AA2=[x*2 for x in AA]
    A1count=map(add, AA2, AB)
    A1F=sum(A1count)/(float(len(AA2))*2)
    return A1F

# output each PRS for each sample, in the form of [sample, *scores],
# and a list of pvalues that are in the order of the scores from each p-value

def labelPRS(PRSdict, sampleNames):
    pList=PRSdict.keys()
    scorelist=[PRSdict[p] for p in pList]
    scoreTable=[sampleScore for sampleScore in zip(sampleNames, *scorelist)]
    return pList, scoreTable

def writePRS(prsTable, outputFile, pvalues, dialect=None):
    title=["Subjects"]
    title.extend(pvalues)
    with open(outputFile, "w") as f:
        csvwriter=csv.writer(f, dialect=dialect)
        csvwriter.writerow(title)
        for score in prsTable:
            csvwriter.writerow(score)
        print("Successfully wrote scores to "+ ntpath.basename(outputFile))




def prepare(gwasfilePath, genofilePath, prunedsnpPath=None):
    gwasfile=sc.textFile(gwasfilePath)
    prunedsnp=sc.textFile(prunedsnpPath)
    mavangeno=sc.textFile(genofilePath)

    gwastable=gwasfile.map(lambda line: line.split("\t"))

    pruneList=prunedsnp.mapPartitions(buildlist).collect()
    pruneSet=set(pruneList)
    pruneLookup=sc.broadcast(pruneSet)

    GWASpruned=gwastable.filter(lambda line: line[0] in pruneLookup.value)  #292232 lines, distinct
    #GWASprunedsnpID=GWASpruned.map(lambda snp: snp[0]).collect()
    #GWASprunedsnpset=set(GWASprunedsnpID)

    GWASprunedsnpTop=GWASpruned.map(lambda snp: (snp[0], snp[3])).collectAsMap()
    print(GWASpruned.take(10))

    genoPruned=mavangeno.filter(lambda line : line.split(' ')[1] in pruneLookup.value).map(lambda x: x.split(" "))
    GWASTopStrandBC=sc.broadcast(GWASprunedsnpTop)
    bpMap=sc.broadcast({"C":"G", "G":"C", "A":"T", "T":"A"})
    MAVANgenotype=genoPruned.map(lambda line: (line[1], makeGenotype(line, GWASTopStrandBC.value, bpMap.value)))
    try:
        test=MAVANgenotype.first()

    except:
        print("Calculation Unccessful")

    return MAVANgenotype, GWASpruned
