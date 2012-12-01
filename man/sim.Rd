\name{sim}

\alias{sim}

\title{Similarity measures for IBD segments}

\description{

  \code{sim}: \R implementation of \code{sim}.

  Similarity measure for IBD segments, tagSNVs, and individuals.

}
\usage{

sim(x,y,simv="minD",minInter=2)

}
\arguments{
  \item{x}{first vector.}
  \item{y}{second vector.}
  \item{simv}{similarity measure: \code{"minD"} (percentage of smaller explained by larger set),
   \code{"jaccard"} (Jaccard index), \code{"dice"} (Dice index), or
   \code{"maxD"}; default \code{"minD"}.}
  \item{minInter}{minimal size of intersection that is required for a non-zero measure; default 2.}

}
\details{

  Similarity measure for IBD segments, tagSNVs, and individuals.

  Implementation in \R.

}
\value{
  \item{erg}{the similarity measure between 0 and 1, where 0 is not
    similar and 1 is identical.}

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

x <- sample(1:15,8)
y <- sample(1:15,8)
sim(x,y,simv="minD",minInter=1)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
