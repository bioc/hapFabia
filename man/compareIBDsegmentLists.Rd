\name{compareIBDsegmentLists}

\alias{compareIBDsegmentLists}
\alias{compareIBDsegmentLists,IBDsegmentList-method}
\alias{compareIBDsegmentLists,IBDsegmentList,ANY,character,ANY,ANY,numeric,numeric-method}

\title{Hierarchical clustering of
  IBD segments stored in IBD segment list(s)}

\description{

  \code{compareIBDsegmentLists}: \R implementation of \code{compareIBDsegmentLists}.

  The IBD segments in one or two list(s) are compared by hierarchical clustering.
  Different similarity measures are available.
  Called by \code{hapFabia}.


}
\usage{

\S4method{compareIBDsegmentLists}{IBDsegmentList,ANY,character,ANY,ANY,numeric,numeric}(IBDsegmentList1,IBDsegmentList2=NULL,simv="minD",pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=6,minIndivid=2)

}

\arguments{
  \item{IBDsegmentList1}{list of IBD segments  given as
    \code{IBDsegmentList} object.}
  \item{IBDsegmentList2}{optional: second list of IBD segments
    given as  \code{IBDsegmentList} object.}
  \item{simv}{similarity measure:
    \code{minD} (percentage of the smaller set explained by the larger set),
    \code{jaccard} (Jaccard index),
    \code{dice} (Dice index), or
    \code{maxD} (percentage of the larger set explained by the smaller set);
    default \code{minD}.}
  \item{pTagSNVs}{optional: exponent for tagSNV similarity.}
  \item{pIndivid}{optional: exponent for individuals similarity.}
  \item{minTagSNVs}{required minimal number of overlapping SNVs to assign
  similarity different from zero; default 6.} 
\item{minIndivid}{required minimal number of
  overlapping individuals to assign
  similarity different from zero; default 2.}
}

\details{

  Similarities are separately computed for SNVs and for individuals using
  one of the similarity measure: \code{"minD"} (percentage of the
  smaller set explained by the larger set), 
  \code{"jaccard"} (Jaccard index), \code{"dice"} (Dice index),
  or  \code{"maxD"}
  (percentage of the larger set explained by the smaller set).
  One minus the product between SNV similarity and
  individuals similarity is the final value
  used for clustering.

  The final similarity measure is not a distance but is symmetric,
  one for similarity of an IBD segment with itself,
  zero if either CNVs or individuals have no overlap,
  and between zero and one.
  Called by \code{hapFabia}.

  Implementation in \R.

}
\value{
  \item{clust}{object of class \code{hclust} which describes the tree
    produced by the hierarchical clustering method.}

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

show(comp)



}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
