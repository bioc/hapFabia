\name{mergeIBDsegmentLists}

\alias{mergeIBDsegmentLists}
\alias{mergeIBDsegmentLists,IBDsegmentList-method}
\alias{mergeIBDsegmentLists,IBDsegmentList,ANY,vector-method}

\title{Merging IBD segments}

\description{

  \code{mergeIBDsegmentLists}: \R implementation of \code{mergeIBDsegmentLists}.

  Merges and combines the IBD segments of one or two list(s).
  A vector gives for each cluster in the
   IBD segment list its new cluster membership as an integer
   number.
  Called by \code{hapFabia} where new membership is determined by
  hierarchical clustering.

}
\usage{

\S4method{mergeIBDsegmentLists}{IBDsegmentList,ANY,vector}(IBDsegmentList1,IBDsegmentList2=NULL,clustIBDsegmentList)


}
\arguments{
  \item{IBDsegmentList1}{object of class \code{IBDsegmentList}.}
  \item{IBDsegmentList2}{optional: second object of class IBDsegmentList.}
 \item{clustIBDsegmentList}{vector giving for each cluster in the
   IBD segment list its new cluster membership as an integer number.}

}
\details{
  A vector gives for each IBD segment its new membership
  as an integer number.
  IBD segments that belong to the same new cluster are merged.

  Implementation in \R.
}
\value{
  \item{IBDsegmentListmerge}{object of \code{IBDsegmentList}
    containing the merged IBD segments.}

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
IBDsegmentList1 <- hapRes$IBDsegmentList1
IBDsegmentList2 <- hapRes$IBDsegmentList2
comp <-
compareIBDsegmentLists(IBDsegmentList1,
   IBDsegmentList2,simv="minD",pTagSNVs=NULL,
   pIndivid=NULL,minTagSNVs=6,minIndivid=2)
if (!is.null(comp)) {
  clustIBDsegmentList <- cutree(comp,h=0.8)
  mergedIBDsegmentList <-
     mergeIBDsegmentLists(IBDsegmentList1=
     IBDsegmentList1,IBDsegmentList2=
     IBDsegmentList2,clustIBDsegmentList=
     clustIBDsegmentList)
}
summary(IBDsegmentList1)
summary(IBDsegmentList2)

summary(mergedIBDsegmentList)
print(IBDsegmentPos(mergedIBDsegmentList[[1]]))
print(IBDsegmentLength(mergedIBDsegmentList[[1]]))

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
