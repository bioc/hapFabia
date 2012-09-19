
#####load libraries#######
library(hapFabia)
library(fabia)

#####convert from .vcf to _mat.txt: step 2. above#######
vcftoFABIA(fileName=fileName)

#####split/ generate segments: step 3. above#######
split_sparse_matrix(fileName=fileName,segmentSize=segmentSize,shiftSize=shiftSize,annotation=TRUE)

#####compute how many segments we have#######
ina <- as.numeric(readLines(paste(fileName,"_mat.txt",sep=""),n=2))
noSNVs <- ina[2]
over <- segmentSize%/%shiftSize
N1 <- noSNVs%/%shiftSize
endRunA <- (N1-over+2)

#####analyze each segment#######
#####may be done by parallel runs#######
iterateSegments(startRun=1,endRun=endRunA,shift=shiftSize,segmentSize=segmentSize,fileName=fileName,individuals=0,upperBP=0.05,p=10,iter=40,alpha=0.03,cyc=50,haploClusterLength=50,Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.03,haplotypes=TRUE)

#####identify duplicates#######
identifyDuplicates(fileName=fileName,startRun=1,endRun=endRunA,shift=shiftSize,segmentSize=segmentSize)

#####analyze results; parallel#######
anaRes <- analyzeHaploClusters(fileName=fileName,startRun=1,endRun=endRunA,shift=shiftSize,segmentSize=segmentSize)
print("Number haplotype clusters:")
print(anaRes$nohaploClusters)
print("Statistics on haplotye cluster length in SNVs (all SNVs in the haplotype cluster):")
print(anaRes$avhaploClusterLengthSNVS)
print("Statistics on haplotye cluster length in bp:")
print(anaRes$avhaploClusterLengthS)
print("Statistics on number of individuals belonging to haplotye clusters:")
print(anaRes$avnoIndividS)
print("Statistics on number of tagSNVs of haplotye clusters:")
print(anaRes$avnoTagSNVsS)
print("Statistics on MAF of tagSNVs of haplotye clusters:")
print(anaRes$avnoFreqS)
print("Statistics on MAF within the group of tagSNVs of haplotye clusters:")
print(anaRes$avnoGroupFreqS)
print("Statistics on number of changes between major and minor allele frequency:")
print(anaRes$avnotagSNVChangeS)
print("Statistics on number of tagSNVs per individual of a haplotype cluster:")
print(anaRes$avnotagSNVsPerIndividualS)
print("Statistics on number of individuals that have the minor allele of tagSNVs:")
print(anaRes$avnoindividualPerTagSNVS)

#####load result for segment 5#######
posAll <- 5 # (5-1)*5000 = 20000: segment 20000 to 21000
start <- (posAll-1)*shiftSize
end <- start + segmentSize
pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
mergedHaploClusterList <- resHapFabia$mergedHaploClusterList # $

#####plot the first specific haplotype cluster in segment 5#######
haploCluster <- mergedHaploClusterList[[1]]
plot(haploCluster,filename=paste(fileName,pRange,"_mat",sep=""))
    ##attention: filename without type ".txt"
