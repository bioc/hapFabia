\name{simu}
\docType{data}
\alias{simu}
\title{Example simulation data for \code{hapFabia}}
\description{

  The data obtained by \code{simulateForFabia} in binary format.
  
}
\usage{simu}
\format{contains a list \code{simu} with  variables
  \enumerate{
    \item \code{namesL},
    \item \code{haploN},
    \item \code{snvs},
    \item \code{annot},
    \item \code{alleleIimp} which contains the information on the
    implanted haplotype cluster: sample names,
  number of chromosomes/individuals,
  number of tagSNVs, annotation of tagSNVs, the genotype data
  where 0 = reference an 1 = minor allele.
}
  \code{annot} is a list with entries:
  \enumerate{
   \item \code{chromosome}, 
   \item \code{position},
   \item \code{snvNames},  
   \item \code{snvMajor},  
   \item \code{snvMinor},  
   \item \code{quality},   
   \item \code{pass},     
   \item \code{info},      
   \item \code{fields},    
   \item \code{frequency}, and 
   \item \code{changed}.
 }
  
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


