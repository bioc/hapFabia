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
    implanted IBD segment: sample names,
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
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datasets}


