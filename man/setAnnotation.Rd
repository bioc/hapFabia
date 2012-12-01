\name{setAnnotation}

\alias{setAnnotation}
\alias{setAnnotation,IBDsegmentList-method}
\alias{setAnnotation,IBDsegmentList,character-method}


\title{Fills in annotations of tagSNVs of a
  list of IBD segments}

\description{
  \code{setAnnotation}: \R implementation of \code{setAnnotation}.

  Fills in the tagSNV annotation of IBD segments
  given in an object of the class \code{IBDsegmentList}.
  The annotation must be given in a file where the first column contains
  the position of the SNV and the second the chromosome.
  The other columns give the
  annotation like "stop gain", "stop loss", "synonymous",
  "non-synonymous", "exonic", "intronic", "intergenic", "promotor",
  etc.
  However other annotations like whether the minor allele is identical
  to the Denisova or Neandertal base can be included.
}

\usage{

\S4method{setAnnotation}{IBDsegmentList,character}(IBDsegmentList,filename)


}
\arguments{
  \item{IBDsegmentList}{object of class \code{IBDsegmentList}.}
  \item{filename}{File containing the SNV annotations, where the first
    column contains the SNV position and the second the chromosome.}


}
\details{

  Implementation in \R.

}
\value{
object of class \code{IBDsegmentList} in which the annotation for the
tagSNVs is set.

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
old_dir <- getwd()
setwd(tempdir())

data(hapRes)

res <- hapRes$res
sPF <- hapRes$sPF
annot <- hapRes$annot
nnL <- length(Z(res)[1,])
labelsA <- cbind(as.character(1:nnL),
   as.character(1:nnL),as.character(1:nnL),
   as.character(1:nnL))
resIBDsegmentList <-
  extractIBDsegments(res=res,sPF=sPF,
  annot=annot,chrom="1",labelsA=labelsA,
  ps=0.9,psZ=0.8,inteA=50,thresA=6,mintagSNVs=6,
  off=0,procMinIndivids=0.1,thresPrune=1e-3)


tagSNVPositions <-
   tagSNVPositions(resIBDsegmentList[[1]])
snvR <- sample(min(tagSNVPositions):max(tagSNVPositions),
   length(tagSNVPositions))
snvA <- sort(unique(c(tagSNVPositions,snvR)))

func = c("stopGain","stopLoss","synonymous",
   "non-synonymous","-","-","-","-","-","-")
for (i in 1:length(snvA)) {

if (i>1) {
write(paste(snvA[i],"1",sample(func,1),sep=" "),
   file="snvAnnotation.txt",ncolumns=100,append=TRUE)
} else {

write(paste(snvA[i],"1",sample(func,1),sep=" "),
   file="snvAnnotation.txt",ncolumns=100,append=FALSE)
}


}

tagSNVAnno(resIBDsegmentList[[1]])

resIBDsegmentList <- setAnnotation(resIBDsegmentList,
   filename="snvAnnotation.txt")

tagSNVAnno(resIBDsegmentList[[1]])

setwd(old_dir)

}


\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
