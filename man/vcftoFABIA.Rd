\name{vcftoFABIA}
\alias{vcftoFABIA}
\title{Converting genotyping data from \code{vcf} to sparse matrix format}
\description{


  \code{vcftoFABIA}: C implementation with an \R wrapper of \code{vcftoFABIA}.

  Converts files giving the genotype in vcf format
  to sparse matrix formats as used by FABIA.
  Phased and unphased genotypes as well as dosages or likelihoods
  are written to files in sparse
  matrix formats.
}
\usage{

vcftoFABIA(fileName,prefixPath="",noSnvs=NULL)

}
\arguments{
  \item{fileName}{string giving the file name without the file type '.vcf'. Attention: remove file type!}
  \item{prefixPath}{path to the file.}
  \item{noSnvs}{optional: the number of SNVs; needed for memory allocation; if not known it is determined by first reading all lines of the file.}

}
\details{


  Implementation in C. Also command line programs are supplied.



}
\seealso{
\code{\link{HaploCluster-class}},
\code{\link{HaploClusterList-class}},
\code{\link{analyzeHaploClusters}},
\code{\link{compareHaploClusterLists}},
\code{\link{extractHaploClusters}},
\code{\link{findDenseRegions}},
\code{\link{hapFabia}},
\code{\link{hapFabiaVersion}},
\code{\link{hapRes}},
\code{\link{chr1ASW1000G}},
\code{\link{haploClusterList2excel}},
\code{\link{identifyDuplicates}},
\code{\link{iterateSegments}},
\code{\link{makePipelineFile}},
\code{\link{matrixPlot}},
\code{\link{mergeHaploClusterLists}},
\code{\link{mergedHaploClusterList}},
\code{\link{plotHaplotypeCluster}},
\code{\link{res}},
\code{\link{setAnnotation}},
\code{\link{setStatistics}},
\code{\link{sim}},
\code{\link{simu}},
\code{\link{simulateHaploClustersFabia}},
\code{\link{simulateHaplotypeClusters}},
\code{\link{split_sparse_matrix}},
\code{\link{toolsFactorizationClass}},
\code{\link{vcftoFABIA}}
}
\author{Sepp Hochreiter}
\examples{


\dontrun{
#########################################
## Already run in "iterateSegments.Rd" ##
#########################################

#Work in a temporary directory.

old_dir <- getwd()
setwd(tempdir())


# Load data and write to vcf file.
data(chr1ASW1000G)
write(chr1ASW1000G,file="chr1ASW1000G.vcf")

#Create the analysis pipeline.
makePipelineFile(fileName="chr1ASW1000G",shiftSize=500,segmentSize=1000)

source("pipeline.R")

# Following files are produced:
list.files(pattern="chr1")



# Next we load segment 5 and there the first and second haplotype cluster
posAll <- 5
start <- (posAll-1)*shiftSize
end <- start + segmentSize
pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
haploClusterList <- resHapFabia$mergedhaploClusterList
summary(haploClusterList)
haploCluster1 <- haploClusterList[[1]]
summary(haploCluster1)
haploCluster2 <- haploClusterList[[2]]
summary(haploCluster2)




#Plot the first haplotype cluster in segment 5
plot(haploCluster1,filename=paste(fileName,pRange,"_mat",sep=""))


#Plot the second haplotype cluster in segment 5
plot(haploCluster2,filename=paste(fileName,pRange,"_mat",sep=""))

setwd(old_dir)

}

\dontrun{
###here an example of the the automatically generated pipeline
### with: shiftSize=5000,segmentSize=10000,fileName="filename"

#####define segments, overlap, filename #######
shiftSize <- 5000
segmentSize <- 10000
fileName="filename" # without type

#####load libraries#######
library(hapFabia)
library(fabia)

#####convert from .vcf to _mat.txt: step 2. above#######
vcftoFABIA(fileName=fileName)

#####split/ generate segments: step 3. above#######
split_sparse_matrix(fileName=fileName,segmentSize=segmentSize,
shiftSize=shiftSize,annotation=TRUE)

#####compute how many segments we have#######
ina <- as.numeric(readLines(paste(fileName,"_mat.txt",sep=""),n=2))
noSNVs <- ina[2]
over <- segmentSize%/%shiftSize
N1 <- noSNVs%/%shiftSize
endRunA <- (N1-over+2)

#####analyze each segment#######
#####may be done by parallel runs#######
iterateSegments(startRun=1,endRun=endRunA,shift=shiftSize,
segmentSize=segmentSize,fileName=fileName,individuals=0,
upperBP=0.05,p=10,iter=40,alpha=0.03,cyc=50,haploClusterLength=50,
Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,
pMAF=0.03,haplotypes=TRUE)

#####identify duplicates#######
identifyDuplicates(fileName=fileName,startRun=1,endRun=endRunA,
shift=shiftSize,segmentSize=segmentSize)

#####analyze results; parallel#######
anaRes <- analyzeHaploClusters(fileName=fileName,startRun=1,endRun=endRunA,
shift=shiftSize,segmentSize=segmentSize)
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

#####load result for segment 50#######
posAll <- 50 # (50-1)*5000 = 245000: segment 245000 to 255000
start <- (posAll-1)*shiftSize
end <- start + segmentSize
pRange <- paste("_",format(start,scientific=FALSE),"_",
format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
haploClusterList <- resHapFabia$mergedHaploClusterList # $

summary(haploClusterList)
#####plot haplotype clusters in segment 50#######
plot(haploClusterList,filename=paste(fileName,pRange,"_mat",sep=""))
   ##attention: filename without type ".txt"

#####plot the first haplotype cluster in segment 50#######

haploCluster <- haploClusterList[[1]]
plot(haploCluster,filename=paste(fileName,pRange,"_mat",sep=""))
   ##attention: filename without type ".txt"

}

}
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{file,connection}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation}
