\name{makePipelineFile}
\alias{makePipelineFile}
\title{makePipelineFile: Make \code{pipleline.R} for  \code{source("Pipleline.R")}}
\description{

  \code{makePipelineFile} creates \code{pipleline.R} for sourcing with
  \code{source("pipleline.R")} to run a whole haplotype cluster extraction pipeline.


}
\usage{

makePipelineFile(fileName,shiftSize=5000,segmentSize=10000) 

}
\arguments{
 \item{fileName}{file name of the genotype file in \code{vcf} format}
 \item{shiftSize}{distance between start of adjacent segments.}
 \item{segmentSize}{number of SNVs or SNPs in a segment.}

}
\details{


  \code{makePipelineFile} creates \code{Pipleline.R} for sourcing with
  \code{source("pipleline.R")} to run a whole haplotype cluster extraction pipeline.

   Attention: this code may run a while for large data sets.

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

old_dir <- getwd()
setwd(tempdir())

makePipelineFile(fileName="genotypeData",
   shiftSize=500,segmentSize=1000)
a <- scan(file = "pipeline.R",
   what = "character")
cat(a)

setwd(old_dir)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
