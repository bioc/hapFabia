\name{split_sparse_matrix}

\alias{split_sparse_matrix}

\title{Splits genotyping data in sparse matrix format into intervals}

\description{

  \code{split_sparse_matrix}:  C implementation with an \R wrapper of \code{split_sparse_matrix}.

  Genotype data of a chromosome is split into intervals where each interval
  leads to a genotype file in sparse matrix format.

}
\usage{

split_sparse_matrix(fileName,sparseMatrixPostfix="_mat.txt",
intervalSize=10000,shiftSize=5000,annotation=TRUE)
}
\arguments{
  \item{fileName}{string giving the genotype data in sparse matrix format
   without type. Attention: no type!}
  \item{sparseMatrixPostfix}{postfix string for sparse matrix format.}
  \item{intervalSize}{number of SNVs in one interval.}
  \item{shiftSize}{number of SNVs between beginning of adjacent intervals
  that is the number of SNVs the intervals are shifted.}
  \item{annotation}{boolean variable indicating whether a annotation file is available.}

}
\details{

  Genotype data in split into intervals of size \code{intervalSize},
  where the distance of the start of adjacent intervals is \code{shiftSize}.
  Thus, it is possible to generate overlapping intervals to account for
  IBD segments that are located at the border of an interval.

  Implementation in C. Also a command line program is supplied.


}
\seealso{
\code{\link{IBDsegment-class}},
\code{\link{IBDsegmentList-class}},
\code{\link{analyzeIBDsegments}},
\code{\link{compareIBDsegmentLists}},
\code{\link{extractIBDsegments}},
\code{\link{findDenseRegions}},
\code{\link{hapFabia}},
\code{\link{hapFabiaVersion}},
\code{\link{hapRes}},
\code{\link{chr1ASW1000G}},
\code{\link{IBDsegmentList2excel}},
\code{\link{identifyDuplicates}},
\code{\link{iterateIntervals}},
\code{\link{makePipelineFile}},
\code{\link{matrixPlot}},
\code{\link{mergeIBDsegmentLists}},
\code{\link{mergedIBDsegmentList}},
\code{\link{plotIBDsegment}},
\code{\link{res}},
\code{\link{setAnnotation}},
\code{\link{setStatistics}},
\code{\link{sim}},
\code{\link{simu}},
\code{\link{simulateIBDsegmentsFabia}},
\code{\link{simulateIBDsegments}},
\code{\link{split_sparse_matrix}},
\code{\link{toolsFactorizationClass}},
\code{\link{vcftoFABIA}}
}
\author{Sepp Hochreiter}
\examples{


\dontrun{
#########################################
## Already run in "iterateIntervals.Rd" ##
#########################################

#Work in a temporary directory.

old_dir <- getwd()
setwd(tempdir())


# Load data and write to vcf file.
data(chr1ASW1000G)
write(chr1ASW1000G,file="chr1ASW1000G.vcf")

#Create the analysis pipeline for IBD segment extraction
makePipelineFile(fileName="chr1ASW1000G",shiftSize=500,intervalSize=1000,haplotypes=TRUE)

source("pipeline.R")

# Following files are produced:
list.files(pattern="chr1")



# Next we load interval 5 and there the first and second IBD segment
posAll <- 5
start <- (posAll-1)*shiftSize
end <- start + intervalSize
pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
IBDsegmentList <- resHapFabia$mergedIBDsegmentList
summary(IBDsegmentList)
IBDsegment1 <- IBDsegmentList[[1]]
summary(IBDsegment1)
IBDsegment2 <- IBDsegmentList[[2]]
summary(IBDsegment2)




#Plot the first IBD segment in interval 5
plot(IBDsegment1,filename=paste(fileName,pRange,"_mat",sep=""))


#Plot the second IBD segment in interval 5
plot(IBDsegment2,filename=paste(fileName,pRange,"_mat",sep=""))

setwd(old_dir)

}

\dontrun{
###here an example of the the automatically generated pipeline
### with: shiftSize=5000,intervalSize=10000,fileName="filename"

#####define intervals, overlap, filename #######
shiftSize <- 5000
intervalSize <- 10000
fileName="filename" # without type
haplotypes <- TRUE
dosage <- FALSE

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

#####split/ generate intervals#######
split_sparse_matrix(fileName=fileName,intervalSize=intervalSize,
shiftSize=shiftSize,annotation=TRUE)

#####compute how many intervals we have#######
ina <- as.numeric(readLines(paste(fileName,"_mat.txt",sep=""),n=2))
noSNVs <- ina[2]
over <- intervalSize\%/\%shiftSize
N1 <- noSNVs\%/\%shiftSize
endRunA <- (N1-over+2)

#####analyze each interval#######
#####may be done by parallel runs#######
iterateIntervals(startRun=1,endRun=endRunA,shift=shiftSize,
intervalSize=intervalSize,fileName=fileName,individuals=0,
upperBP=0.05,p=10,iter=40,alpha=0.03,cyc=50,IBDsegmentLength=50,
Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,
pMAF=0.03,haplotypes=haplotypes,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,
simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)

#####identify duplicates#######
identifyDuplicates(fileName=fileName,startRun=1,endRun=endRunA,
shift=shiftSize,intervalSize=intervalSize)

#####analyze results; parallel#######
anaRes <- analyzeIBDsegments(fileName=fileName,startRun=1,endRun=endRunA,
shift=shiftSize,intervalSize=intervalSize)
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

#####load result for interval 50#######
posAll <- 50 # (50-1)*5000 = 245000: interval 245000 to 255000
start <- (posAll-1)*shiftSize
end <- start + intervalSize
pRange <- paste("_",format(start,scientific=FALSE),"_",
format(end,scientific=FALSE),sep="")
load(file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))
IBDsegmentList <- resHapFabia$mergedIBDsegmentList # $

summary(IBDsegmentList)
#####plot IBD segments in interval 50#######
plot(IBDsegmentList,filename=paste(fileName,pRange,"_mat",sep=""))
   ##attention: filename without type ".txt"

#####plot the first IBD segment in interval 50#######

IBDsegment <- IBDsegmentList[[1]]
plot(IBDsegment,filename=paste(fileName,pRange,"_mat",sep=""))
   ##attention: filename without type ".txt"

}

}
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\value{Splits genotyping data in sparse matrix format into intervals}
\keyword{file,connection}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation}
