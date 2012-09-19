\name{compareHaploClusterLists}

\alias{compareHaploClusterLists}
\alias{compareHaploClusterLists,HaploClusterList-method}
\alias{compareHaploClusterLists,HaploClusterList,ANY,character,ANY,ANY,numeric,numeric-method}

\title{Hierarchical clustering of haplotype
  clusters stored in haplotype cluster list(s)}

\description{

  \code{compareHaploClusterLists}: \R implementation of \code{compareHaploClusterLists}.

  The haplotype clusters in one or two list(s) are compared by hierarchical clustering.
  Different similarity measures are available.
  Called by \code{hapFabia}.


}
\usage{

\S4method{compareHaploClusterLists}{HaploClusterList,ANY,character,ANY,ANY,numeric,numeric}(haploClusterList1,haploClusterList2=NULL,simv="minD",pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=6,minIndivid=2)

}

\arguments{
  \item{haploClusterList1}{list of haplotype clusters  given as
    \code{HaploClusterList} object.}
  \item{haploClusterList2}{optional: second list of haplotype clusters
    given as  \code{HaploClusterList} object.}
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
  one for similarity of a haplotype cluster with itself,
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

show(comp)



}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
