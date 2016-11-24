# convert form .gen format in to single digit format format
from operator import add
def makeGenotype(line,idCol):
    AA=map(float,line[5::3])
    AB=map(float, line[6::3])
    AA2=[x*2 for x in AA]
    genotype=[AA2[i]+AB[i] for i in range(len(AA2))]
    return (line[idCol], genotype)

bpPair={"C":"G", "G":"C", "A":"T", "T":"A"}

def makeGenotypeCheckRef(line, freqMap, gwasMap, idCol, a1Col, startCol, gwasA1, gwasMAF, bpMap=bpPair):
    rsid=line[idCol]
    genosnp = (line[a1Col], line[a1Col+1])
    gwassnp = (gwasMap[rsid][gwasA1],gwasMap[rsid][gwasA1+1])

    if rsid in gwasMap:
        if genosnp[0] != bpMap[genosnp[1]]:   # if the SNP is not ambiguous, a.k.a not (A,T) or (C,G) pair
            if genosnp[0] == gwassnp[0] or bpMap[genosnp[0]]==gwassnp[0]:
                AA=list(map(float,line[startCol::3]))
                AB=list(map(float, line[(startCol+1)::3]))
                AA2=[x*2 for x in AA]
                genotype=list(map(add, AA2, AB))
                return ((rsid, line[a1Col], line[a1Col+1]), genotype)

            elif genosnp[1] == gwassnp[0] or bpMap[genosnp[1]]==gwassnp[0]:
                AA=list(map(float, line[(startCol+2)::3]))
                AB=list(map(float, line[(startCol+1)::3]))
                AA2=[x*2 for x in AA]
                genotype=list(map(add, AA2, AB))
                return ((rsid, line[a1Col], line[a1Col+1]), genotype)
            else: pass

        elif genosnp[0] == bpMap[genosnp[1]]:  # if the snp is ambiguous
            genofreq=freqMap[rsid]
            majorIdx=genofreq.index(max(genofreq))
            gwasMajorIdx=1-int(gwasMap[rsid][gwasMAF]>0.5)
            if majorIdx==gwasMajorIdx:
                AA=list(map(float, line[startCol::3]))
            else:
                AA=list(map(float, line[(startCol+2)::3]))

            AB=list(map(float, line[6::3]))
            AA2=[x*2 for x in AA]
            genotype=list(map(add, AA2, AB))

            return ((rsid, line[a1Col], line[a1Col+1]), genotype)

    else:
        pass



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

# calculate PRS from genotype
def calcPRSFromGeno(genotypeRDD, oddsMap):
    oddsBC=sc.broadcast(oddsMap)
    totalcount=genotypeRDD.count()
    multiplied=genotypeRDD.map(lambda line:multiplyOdds(line, oddsMap))
    filtered=multiplied.filter(lambda line: line is not None)
    PRS=multiplied.reduce(lambda a,b: [a[i]+b[i] for i in range(len(a))])
    normalizedPRS=[x/totalcount for x in PRS]
    return PRS

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

def getMaf(geno, rsidCol, start, a1):
    AA=line[0::3]
    AB=line[1::3]
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
