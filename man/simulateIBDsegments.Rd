\name{simulateIBDsegments}

\alias{simulateIBDsegments}

\title{Generates simulated genotyping data with IBD segments}

\description{

  \code{simulateIBDsegments}: \R implementation of \code{simulateIBDsegments}.

  Genotype data with rare variants is simulated.
  Into these datan IBD segments are implanted.
   All data sets and information are written to files.
}
\usage{

simulateIBDsegments(fileprefix="dataSim",minruns=1,
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
  \item{noImplanted}{number of IBD segments that are implanted.}
  \item{implanted}{number of individuals belonging to specific IBD segment.}
  \item{length}{length of the IBD segments in number of SNVs.}
  \item{minors}{number of tagSNVs for each IBD segment.}
  \item{mismatches}{number of minor allele tagSNV mismatches
    for individuals belonging
    to the IBD segment.}
  \item{mismatchImplanted}{percentage of individuals of
    an IBD segment that have mismatches.}
  \item{overlap}{minimal overlap of the founder interval between
    individuals belonging to a specific IBD segment (the interval
    may be broken at the ends).}
  \item{noOverwrite}{\code{noOverwrite=TRUE} ensures that an IBD segment is not superimposed by another IBD segment.}


}
\details{

  Data simulations focuses on rare variants but common variants are
  possible, too.
  Linkage disequilibrium and haplotype blocks are not simulated except by implanting IBD segments.

  Simulated data is written to files. For BEAGLE the data is written to "...beagle.txt".
  For PLINK the data is written to "...plink.ped", "...plink.map", and "...plink.fam".
  For the MCMC method the data is written to "...mcmc.genotype", "...mcmc.posmaf", and "...mcmc.initz".
  For RELATE the data is written to "...relate.geno", "...relate.pos", and "...relate.chr".
  For \pkg{fabia} the data is written
  to "...fabia_individuals.txt", "...fabia_annot.txt" "...fabia_mat.txt".

   Information on parameters for data simulation is written to
   "...Parameters.txt" while information on implanted IBD segments is written to "...Impl.txt".

   Most information is also written in \R binary ".Rda" files.

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
\dontrun{
old_dir <- getwd()
setwd(tempdir())

simulateIBDsegments(minruns=1,maxruns=1,snvs=1000,individualsN=10,avDistSnvs=100,avDistMinor=15,noImplanted=1,implanted=10,length=100,minors=10,mismatches=0,mismatchImplanted=0.5,overlap=50,noOverwrite=FALSE) 

setwd(old_dir)

}

}
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datagen}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
