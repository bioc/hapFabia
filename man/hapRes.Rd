\name{hapRes}
\docType{data}
\alias{hapRes}
\title{hapFabia Data: hapRes (result of \code{hapFabia}) }
\description{

  Results of a \code{hapFabia} call.
  List containing
 \enumerate{
 \item \code{mergedhaploClusterList}: an object of the class
   \code{HaploClusterList} that contains the extracted haplotype clusters
   that were extracted from two histograms with different offset.
 \item \code{res}: the result of FABIA.
 \item \code{sPF}: samples per loading of this FABIA result.
 \item \code{annot}: annotation for the genotype data.
 \item \code{haploClusterList1}: an object of the class
   \code{HaploClusterList} that contains the result of haplotype cluster extraction from the first histogram.
 \item \code{haploClusterList2}: an object of the class
   \code{HaploClusterList} that contains the result of haplotype cluster extraction from the second histogram.
 \item \code{mergedhaploClusterList1}: an object of the class
   \code{HaploClusterList} that contains the merged result of the first haplotype cluster
   extraction (redundancies removed).
 \item \code{mergedhaploClusterList2}: an object of the class
   \code{HaploClusterList} that contains the merged result of the second haplotype cluster
   extraction (redundancies removed).

}
  
}
\usage{hapRes}
\format{  List containing
 \enumerate{
 \item \code{mergedhaploClusterList}: an object of the class
   \code{HaploClusterList} that contains the extracted haplotype clusters
   that were extracted from two histograms with different offset.
 \item \code{res}: the result of FABIA.
 \item \code{sPF}: samples per loading of this FABIA result.
 \item \code{annot}: annotation for the genotype data.
 \item \code{haploClusterList1}: an object of the class
   \code{HaploClusterList} that contains the result of haplotype cluster extraction from the first histogram.
 \item \code{haploClusterList2}: an object of the class
   \code{HaploClusterList} that contains the result of haplotype cluster extraction from the second histogram.
 \item \code{mergedhaploClusterList1}: an object of the class
   \code{HaploClusterList} that contains the merged result of the first haplotype cluster
   extraction (redundancies removed).
 \item \code{mergedhaploClusterList2}: an object of the class
   \code{HaploClusterList} that contains the merged result of the second haplotype cluster
   extraction (redundancies removed).
}
}
\source{result of a call of \code{hapFabia} }
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
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datasets}


