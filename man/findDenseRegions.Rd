\name{findDenseRegions}
\alias{findDenseRegions}
\title{findDenseRegions: Find accumulations of values via histogram counts}
\description{

  \code{findDenseRegions}: \R implementation of \code{findDenseRegions}.

  Extracts histogram bins which counts are larger than a threshold.
  Only L values larger than a given quantile are considered.
}
\usage{

findDenseRegions(obs,p=0.9,inte=500,thres=11,off=0)


}
\arguments{
  \item{obs}{values for constructing the histogram.}
 \item{p}{quantile above which the values of \code{obs} are included into the histogram.}
  \item{inte}{size of the histogram bins.}
  \item{thres}{threshold for histogram counts: bin with counts larger equal the threshold are selected.}
  \item{off}{offset of the histogram bin positions.}
}
\details{

  Extracts histogram bins which counts are larger than a threshold.
  The threshold is supplied and can be computed according to some assumptions
  on expected bin counts.
  Only L values larger than a given quantile are considered.



  Implementation in \R.

}
\value{
list with
  \item{m}{vector of locations of the selected bin (middle of bins).}
  \item{l}{vector of lengths of the bins.}
  \item{pos}{list where each element is a vector of locations of values that contributed to the counts (SNV positions).}
  \item{len}{vector of counts or equivalently vector of the number of SNVs.}

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

data(res)
ib <- findDenseRegions(L(res)[,1],p=0.9,
   inte=50,thres=6,off=0)
print(ib$len)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
