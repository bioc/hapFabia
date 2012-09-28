\name{setStatistics}

\alias{setStatistics}
\alias{setStatistics,HaploClusterList-method}


\title{Computes and stores the statistics of a haplotype cluster
  list}

\description{

  \code{setStatistics}: \R implementation of \code{setStatistics}.

  Computes the statistics of a haplotype cluster
  list given as an object of the class \code{HaploClusterList}.
  In the slot \code{statistics} of an object of the class
  \code{HaploClusterList} the summary statistics
  across the list of haplotype clusters are stored.
  Following characteristics are stored in the list \code{statistics}:
  \enumerate{
    \item "avhaploClusterPosS": physical position
\item "avhaploClusterLengthSNVS": length in SNVs
\item "avhaploClusterLengthS": length in bp
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

\S4method{setStatistics}{HaploClusterList}(haploClusterList)


}
\arguments{
  \item{haploClusterList}{object of class \code{HaploClusterList}.}

}
\details{

  The list can be extended by the user.

  Implementation in \R.

}
\value{
object of class  \code{HaploClusterList} with statistics set
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
res <- hapRes$res
sPF <- hapRes$sPF
annot <- hapRes$annot
nnL <- length(Z(res)[1,])
labelsA <- cbind(as.character(1:nnL),as.character(1:nnL),
   as.character(1:nnL),as.character(1:nnL))
resHaploClusterList <-
  extractHaploClusters(res=res,sPF=sPF,annot=annot,
  chrom="1",labelsA=labelsA,ps=0.9,psZ=0.8,inteA=50,
  thresA=6,mintagSNVs=6,off=0,procMinIndivids=0.1,
  thresPrune=1e-3)

summary(resHaploClusterList)

resHaploClusterList <- setStatistics(resHaploClusterList)

summary(resHaploClusterList)


}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}