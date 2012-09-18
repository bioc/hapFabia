\name{simu}
\docType{data}
\alias{simu}
\title{hapFabia Data: simulation (simulation data for \code{hapFabia}) }
\description{

  The data obtaind by \code{simulateForFabia} in binary to
  write to files.
  
}
\usage{simu}
\format{contains a list \code{simu} with  variables \code{namesL}, \code{haploN},
  \code{snvs}, \code{annot}, \code{alleleIimp} which are sample names,
  number of chromosomes/individuals,
  number of tagSNVs, annotation of tagSNVs, the data
  where 0 = reference an 1 = minor allele.

  \code{annot} is a list with entries:
  \code{chromosome}, 
  \code{position},
  \code{snvNames},  
  \code{snvMajor},  
  \code{snvMinor},  
  \code{quality},   
  \code{pass},     
  \code{info},      
  \code{fields},    
  \code{frequency}, and 
  \code{changed}.
   
  }
 
\source{from \code{simulateForFabia}}
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


