\name{simulateIBDsegmentsFabia}

\alias{simulateIBDsegmentsFabia}

\title{Generates simulated genotyping data with IBD segments for \pkg{fabia}}

\description{

  \code{simulateIBDsegmentsFabia}: \R implementation of \code{simulateIBDsegmentsFabia}.

   Genotype data is simulated which contains rare variants and
   implanted IBD segments.
   Output is written for the bicluster algorithm \pkg{fabia}.


}
\usage{

simulateIBDsegmentsFabia(fileprefix="dataSim",
   minruns=1,maxruns=1,snvs=1000,individualsN=100,
   avDistSnvs=100,avDistMinor=10,noImplanted=1,
   implanted=10,length=50,minors=30,mismatches=0,
   mismatchImplanted=0.5,overlap=50)

}
\arguments{
  \item{fileprefix}{prefix of file names containing data generated in this simulation.}
  \item{minruns}{start index for generating multiple data sets.}
  \item{maxruns}{end index for generating multiple data sets.}
  \item{snvs}{number of SNVs in this simulation.}
  \item{individualsN}{number of individuals in this simulation.}
  \item{avDistSnvs}{average genomic distance in bases between SNVs.}
  \item{avDistMinor}{average distance between minor alleles, thus \code{1/avDistMinor}
      is the average minor allele frequency (MAF).}
  \item{noImplanted}{number of IBD segments that are implanted.}
  \item{implanted}{number of individuals into which a specific IBD segment is implanted.}
  \item{length}{length of the IBD segments in number of SNVs.}
  \item{minors}{number of tagSNVs for each IBD segment.}
  \item{mismatches}{number of base mismatches of an implanted IBD segment to the original IBD segment.}
  \item{mismatchImplanted}{percentage of IBD segment occurrence that have mismatches.}
  \item{overlap}{minimal IBD segment overlap between implanted IBD segments (they are broken at the ends).}

}
\details{

  Data simulations for \pkg{fabia} focuses on rare variants but common variants are possible.
  Linkage disequilibrium and haplotype blocks are not simulated except by implanting IBD segments.

  Simulated data is written
  to "...fabia_individuals.txt", "...fabia_annot.txt" "...fabia_mat.txt".


  Implementation in \R.

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

simulateIBDsegmentsFabia()

setwd(old_dir)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datagen}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
