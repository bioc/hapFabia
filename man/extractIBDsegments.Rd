\name{extractIBDsegments}
\alias{extractIBDsegments}
\alias{extractIBDsegments,Factorization-method}
\alias{extractIBDsegments,Factorization,list,data.frame,character,matrix,numeric,numeric,numeric,numeric,numeric,numeric,numeric,numeric-method}



\title{Extract IBD segments from a
  \code{fabia} result}

\description{

  \code{extractIBDsegments}: \R implementation of \code{extractIBDsegments}.

  IBD segments are identified in FABIA \code{Factorization} objects.
  First accumulations of correlated SNVs are found.
  Then IBD segments in these accumulations are disentangled.
  Finally IBD segments are pruned off spurious correlated SNVs.


}
\usage{

\S4method{extractIBDsegments}{Factorization,list,data.frame,character,matrix,numeric,numeric,numeric,numeric,numeric,numeric,numeric,numeric}(res,sPF,annot=NULL,chrom="",labelsA=NULL,ps=0.9,psZ=0.8,inteA=500,thresA=11,mintagSNVs=8,off=0,procMinIndivids=0.1,thresPrune=1e-3)


}
\arguments{
  \item{res}{result of \code{fabia} given as \code{Factorization} object.}
  \item{sPF}{genotype data obtained by \code{fabia} procedure
    \code{samplesPerFeature}; it gives for each SNV the
    individuals/chromosomes that possess the minor allele.} 
  \item{annot}{annotation for the tagSNVs as an object of the class
    \code{data.frame}; if it is \code{NULL} then a dummy annotation is
    generated.} 
  \item{chrom}{the chromosome the genotyping data stems from.}
  \item{labelsA}{labels for the individuals; if it is \code{NULL} then
    dummy labels by enumerating individuals are generated.} 
 \item{ps}{quantile above which the L values are considered for
   IBD segment extraction.} 
  \item{psZ}{quantile above which the largest Z values are considered
    for IBD segment extraction.} 
  \item{inteA}{number of SNVs in a histogram bin which correspond to the
    desired IBD segment length.}
  \item{thresA}{threshold for histogram counts above which SNVs are
    viewed to be locally accumulated in a histogram bin.}
 \item{mintagSNVs}{threshold for minimal tagSNV overlap of intervals in a
   IBD segment.} 
   \item{off}{offset of the histogram.}
   \item{procMinIndivids}{percent of cluster individuals that must have
     the minor allele to consider an SNV as IBD segment tagSNV.}
   \item{thresPrune}{threshold on the probability of having minimal
     distance to neighboring tagSNVs; used to prune off SNVs at the
     border of IBD segments.} 
}
\details{


  The threshold  \code{thresA} for counts in a bin, which indicates
  SNV accumulations, is computed and provided by
  \code{hapFabia} when calling this method.
  Distance probabilities for pruning are based on an exponential
  distribution with 
  the median distance between tagCNVs as parameter (one over the rate).
  Thus, the counts are assumed to be Poisson distributed.
  At the IBD segment border, SNVs that have a large
  distance to the closest tagSNV are pruned off.
  \code{thresPrune} gives the pruning threshold via a \eqn{p}{p}-value for
  observing this distance or a larger based on the exponential distribution.

  Implementation in \R.

}
\value{
An instance of the class \code{IBDsegmentList} containing the
extracted IBD segments.
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
labelsA <- cbind(as.character(1:nnL),
   as.character(1:nnL),as.character(1:nnL),
   as.character(1:nnL))
resIBDsegmentList <- extractIBDsegments(res=res,
   sPF=sPF,annot=annot,chrom="1",labelsA=labelsA,
   ps=0.9,psZ=0.8,inteA=50,thresA=6,mintagSNVs=6,
   off=0,procMinIndivids=0.1,thresPrune=1e-3)

summary(resIBDsegmentList)

print("Position of the first IBD segment:")
print(IBDsegmentPos(resIBDsegmentList[[1]]))

print("Length of the first IBD segment:")
print(IBDsegmentLength(resIBDsegmentList[[1]]))


}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNVs,rare IBD segments,short IBD segments}

