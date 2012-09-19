\name{sim}

\alias{sim}

\title{Similarity measures for haplotype clusters}

\description{

  \code{sim}: \R implementation of \code{sim}.

  Similarity measure for haplotype clusters, tagSNVs, and individuals.

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

  Similarity measure for haplotype clusters, tagSNVs, and individuals.

  Implementation in \R.

}
\value{
  \item{erg}{the similarity measure between 0 and 1, where 0 is not
    similar and 1 is identical.}

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
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
