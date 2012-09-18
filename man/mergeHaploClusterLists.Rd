\name{mergeHaploClusterLists}

\alias{mergeHaploClusterLists}
\alias{mergeHaploClusterLists,HaploClusterList-method}
\alias{mergeHaploClusterLists,HaploClusterList,ANY,vector-method}

\title{mergeHaploClusterLists: Merging and combining haplotype clusters or IBD tracts}

\description{

  \code{mergeHaploClusterLists}: \R implementation of \code{mergeHaploClusterLists}.

  Merges and combines the haplotype clusters of one or two list(s).
  A vector vector gives for each cluster in the
   haplotype cluster list its new cluster membership as an interger
   number.
  Called by \code{hapFabia}.

}
\usage{

\S4method{mergeHaploClusterLists}{HaploClusterList,ANY,vector}(haploClusterList1,haploClusterList2=NULL,clustHaploClustList)


}
\arguments{
  \item{haploClusterList1}{object of class \code{HaploClusterList}.}
  \item{haploClusterList2}{optional: second object of class HaploClusterList.}
 \item{clustHaploClustList}{vector giving for each cluster in the
   haplotype cluster list its new cluster membership as an interger number.}

}
\details{
  A vector givesfor each haplotype clusters its membership
  to a new cluster as an integer number.
  Haplotype clusters that belong to the same new cluster are merged,
  where for example individuals
  and tagSNVs of each cluster are merged.

  Implementation in \R.
}
\value{
  \item{haploClusterListmerge}{object of \code{HaploClusterList}
    containing the merged haplotype clusters.}

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


data(hapRes)
haploClusterList1 <- hapRes$haploClusterList1
haploClusterList2 <- hapRes$haploClusterList2
comp <-
compareHaploClusterLists(haploClusterList1,
   haploClusterList2,simv="minD",pTagSNVs=NULL,
   pIndivid=NULL,minTagSNVs=6,minIndivid=2)
if (!is.null(comp)) {
  clustHaploClustList <- cutree(comp,h=0.8)
  mergedhaploClusterList <-
     mergeHaploClusterLists(haploClusterList1=
     haploClusterList1,haploClusterList2=
     haploClusterList2,clustHaploClustList=
     clustHaploClustList)
}
summary(haploClusterList1)
summary(haploClusterList2)

summary(mergedhaploClusterList)
print(haploClusterPos(mergedhaploClusterList[[1]]))
print(haploClusterLength(mergedhaploClusterList[[1]]))

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
