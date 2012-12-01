
#####load library#######
library(hapFabia)

#####convert from .vcf to _mat.txt#######
vcftoFABIA(fileName=fileName)

#####copy haplotype, genotype, or dosage matrix to matrix#######
if (haplotypes) {
    file.copy(paste(fileName,"_matH.txt",sep=""), paste(fileName,"_mat.txt",sep=""))
} else {
    if (dosage) {
        file.copy(paste(fileName,"_matD.txt",sep=""), paste(fileName,"_mat.txt",sep=""))
    } else {
        file.copy(paste(fileName,"_matG.txt",sep=""), paste(fileName,"_mat.txt",sep=""))
    }
}

#####split/ generate segments#######
split_sparse_matrix(fileName=fileName,intervalSize=intervalSize,shiftSize=shiftSize,annotation=TRUE)

#####compute how many segments we have#######
ina <- as.numeric(readLines(paste(fileName,"_mat.txt",sep=""),n=2))
noSNVs <- ina[2]
over <- intervalSize%/%shiftSize
N1 <- noSNVs%/%shiftSize
endRunA <- (N1-over+2)

#####analyze each segment#######
#####may be done by parallel runs#######
iterateIntervals(startRun=1,endRun=endRunA,shift=shiftSize,intervalSize=intervalSize,fileName=fileName,individuals=0,upperBP=0.05,p=10,iter=40,alpha=0.03,cyc=50,IBDsegmentLength=50,Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.03,haplotypes=haplotypes,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,
simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)

#####identify duplicates#######
identifyDuplicates(fileName=fileName,startRun=1,endRun=endRunA,shift=shiftSize,intervalSize=intervalSize)

#####analyze results; parallel#######
anaRes <- analyzeIBDsegments(fileName=fileName,startRun=1,endRun=endRunA,shift=shiftSize,intervalSize=intervalSize)
print("Number IBD segments:")
print(anaRes$noIBDsegments)
print("Statistics on IBD segment length in SNVs (all SNVs in the IBD segment):")
print(anaRes$avIBDsegmentLengthSNVS)
print("Statistics on IBD segment length in bp:")
print(anaRes$avIBDsegmentLengthS)
print("Statistics on number of individuals belonging to IBD segments:")
print(anaRes$avnoIndividS)
print("Statistics on number of tagSNVs of IBD segments:")
print(anaRes$avnoTagSNVsS)
print("Statistics on MAF of tagSNVs of IBD segments:")
print(anaRes$avnoFreqS)
print("Statistics on MAF within the group of tagSNVs of IBD segments:")
print(anaRes$avnoGroupFreqS)
print("Statistics on number of changes between major and minor allele frequency:")
print(anaRes$avnotagSNVChangeS)
print("Statistics on number of tagSNVs per individual of an IBD segment:")
print(anaRes$avnotagSNVsPerIndividualS)
print("Statistics on number of individuals that have the minor allele of tagSNVs:")
print(anaRes$avnoindividualPerTagSNVS)

#####load result for segment 5#######
posAll <- 5 # (5-1)*5000 = 20000: segment 20000 to 21000
start <- (posAll-1)*shiftSize
end <- start + intervalSize
pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
mergedIBDsegmentList <- resHapFabia$mergedIBDsegmentList # $

#####plot the first specific IBD segment in segment 5#######
IBDsegment <- mergedIBDsegmentList[[1]]
plot(IBDsegment,filename=paste(fileName,pRange,"_mat",sep=""))
    ##attention: filename without type ".txt"
