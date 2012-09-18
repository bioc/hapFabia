\name{simulateHaplotypeClusters}

\alias{simulateHaplotypeClusters}

\title{simulateHaplotypeClusters: Constructs simulation data with rare variants and implanted haplotype clusters}

\description{

  \code{simulateHaplotypeClusters}: \R implementation of \code{simulateHaplotypeClusters}.

  Genotype data with rare variants is simulated.
  Into these data haplotype clusters are implanted.
   All data sets and information are written to files.
}
\usage{

simulateHaplotypeClusters(fileprefix="dataSim",minruns=1,
   maxruns=100,snvs=10000,individualsN=100,avDistSnvs=100,
   avDistMinor=25,noImplanted=1,implanted=10,length=100,
   minors=20,mismatches=0,mismatchImplanted=0.5,overlap=50,
   noOverwrite=FALSE) 


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
  \item{implanted}{number of individuals belonging to specific haplotype cluster.}
  \item{length}{length of the haplotype clusters in number of SNVs.}
  \item{minors}{number of tagSNVs for each haplotype cluster.}
  \item{mismatches}{number of minor allele tagSNV mismatches
    for individuals belonging
    to the haplotype cluster.}
  \item{mismatchImplanted}{percentage of individuals of a haplotype
    cluster that have mismatches.}
  \item{overlap}{minimal overlap of the founder segment between
    individuals belonging to a specific haplotype cluster (the segment
    may be broken at the ends).}
  \item{noOverwrite}{\code{noOverwrite=TRUE} ensures that a haplotype cluster is not superimposed by another haplotype cluster.}


}
\details{

  Data simulations focuses on rare variants but common variants are
  possible, too.
  Linkage disequilibrium and haplotype blocks are not simulated except by implanting haplotype clusters.

  Simulated data is written to files. For BEAGLE the data is written to "...beagle.txt".
  For PLINK the data is written to "...plink.ped", "...plink.map", and "...plink.fam".
  For the MCMC method the data is written to "...mcmc.genotype", "...mcmc.posmaf", and "...mcmc.initz".
  For RELATE the data is written to "...relate.geno", "...relate.pos", and "...relate.chr".
  For \pkg{fabia} the data is written
  to "...fabia_individuals.txt", "...fabia_annot.txt" "...fabia_mat.txt".

   Information on parameters for data simulation is written to
   "...Parameters.txt" while information on implanted haplotype clusters is written to "...Impl.txt".

   Most information is also written in \R binary ".Rda" files.

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
\dontrun{
old_dir <- getwd()
setwd(tempdir())

simulateHaplotypeClusters(minruns=1,maxruns=1,snvs=1000,individualsN=10,avDistSnvs=100,avDistMinor=15,noImplanted=1,implanted=10,length=100,minors=10,mismatches=0,mismatchImplanted=0.5,overlap=50,noOverwrite=FALSE) 

setwd(old_dir)

}

}
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datagen}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
