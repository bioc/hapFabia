\name{setStatistics}

\alias{setStatistics}
\alias{setStatistics,IBDsegmentList-method}


\title{Computes and stores the statistics of an IBD segment
  list}

\description{

  \code{setStatistics}: \R implementation of \code{setStatistics}.

  Computes the statistics of an IBD segment
  list given as an object of the class \code{IBDsegmentList}.
  In the slot \code{statistics} of an object of the class
  \code{IBDsegmentList} the summary statistics
  across the list of IBD segments are stored.
  Following characteristics are stored in the list \code{statistics}:
  \enumerate{
    \item "avIBDsegmentPosS": physical position
\item "avIBDsegmentLengthSNVS": length in SNVs
\item "avIBDsegmentLengthS": length in bp
\item "avnoIndividS": number individuals
\item "avnoTagSNVsS": number tagSNVs
\item "avnoFreqS": tagSNV frequency
\item "avnoGroupFreqS": tagSNV group frequency
\item "avnotagSNVChangeS":  tagSNV change between
  minor and major allele
\item "avnotagSNVsPerIndividualS": tagSNVs per individual
\item "avnoindividualPerTagSNVS": individuals per
  tagSNV.
  }
}

\usage{

\S4method{setStatistics}{IBDsegmentList}(IBDsegmentList)


}
\arguments{
  \item{IBDsegmentList}{object of class \code{IBDsegmentList}.}

}
\details{

  The list can be extended by the user.

  Implementation in \R.

}
\value{
object of class  \code{IBDsegmentList} with statistics set
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
data(hapRes)
res <- hapRes$res
sPF <- hapRes$sPF
annot <- hapRes$annot
nnL <- length(Z(res)[1,])
labelsA <- cbind(as.character(1:nnL),as.character(1:nnL),
   as.character(1:nnL),as.character(1:nnL))
resIBDsegmentList <-
  extractIBDsegments(res=res,sPF=sPF,annot=annot,
  chrom="1",labelsA=labelsA,ps=0.9,psZ=0.8,inteA=50,
  thresA=6,mintagSNVs=6,off=0,procMinIndivids=0.1,
  thresPrune=1e-3)

summary(resIBDsegmentList)

resIBDsegmentList <- setStatistics(resIBDsegmentList)

summary(resIBDsegmentList)


}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
