\name{simulateHaploClustersFabia}

\alias{simulateHaploClustersFabia}

\title{simulateHaploClustersFabia:  Constructs simulation data with rare variants and implanted haplotype clusters for \pkg{fabia}}

\description{

  \code{simulateHaploClustersFabia}: \R implementation of \code{simulateHaploClustersFabia}.

   For bicluster algorithm \pkg{fabia}
   genotype data with rare variants is simulated and
   subsequently haplotype clusters implanted.

}
\usage{

simulateHaploClustersFabia(fileprefix="dataSim",
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
  \item{noImplanted}{number of haplotype clusters that are implanted.}
  \item{implanted}{number of individuals into which a specific haplotype cluster is implanted.}
  \item{length}{length of the haplotype clusters in number of SNVs.}
  \item{minors}{number of tagSNVs for each haplotype cluster.}
  \item{mismatches}{number of base mismatches of an implanted haplotype cluster to the original haplotype cluster.}
  \item{mismatchImplanted}{percentage of haplotype cluster occurrence that have mismatches.}
  \item{overlap}{minimal haplotype cluster overlap between implanted haplotype clusters (they are broken at the ends).}

}
\details{

  Data simulations for \pkg{fabia} focuses on rare variants but common variants are possible.
  Linkage disequilibrium and haplotype blocks are not simulated except by implanting haplotype clusters.

  Simulated data is written
  to "...fabia_individuals.txt", "...fabia_annot.txt" "...fabia_mat.txt".


  Implementation in \R.

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

simulateHaploClustersFabia()

setwd(old_dir)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datagen}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
