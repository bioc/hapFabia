\name{makePipelineFile}
\alias{makePipelineFile}
\title{Generate \code{pipleline.R}}
\description{

  \code{makePipelineFile} creates \code{pipleline.R} for sourcing with
  \code{source("pipleline.R")} to run a whole IBD segment
  extraction pipeline.


}
\usage{

makePipelineFile(fileName,shiftSize=5000,intervalSize=10000,haplotypes=FALSE,dosage=FALSE) 

}
\arguments{
 \item{fileName}{file name of the genotype file in \code{vcf} format}
 \item{shiftSize}{distance between start of adjacent intervals.}
 \item{intervalSize}{number of SNVs in a interval.}
 \item{haplotypes}{should haplotypes (phased genotypes) be used; default
   FALSE.}
 \item{dosage}{should dosages be used if haplotypes is FALSE; default FALSE.}

}
\details{


  \code{makePipelineFile} creates \code{Pipleline.R} for sourcing with
  \code{source("pipleline.R")} to run a whole IBD segment extraction pipeline.

   Attention: this code may run a while for large data sets.

}
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



\author{Sepp Hochreiter}

\examples{

old_dir <- getwd()
setwd(tempdir())

makePipelineFile(fileName="genotypeData",
   shiftSize=500,intervalSize=1000)
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
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
