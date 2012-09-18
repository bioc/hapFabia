\name{plotHaplotypeCluster}

\alias{plotHaplotypeCluster}

\title{plotHaplotypeCluster: Plots a haplotype cluster given individuals, tagSNVs, and their physical positions}

\description{

  \code{plotHaplotypeCluster}: \R implementation of \code{plotHaplotypeCluster}.

  A haplotype cluster is plotted where individuals
  to plot are provided together with tagSNVs and their
  physical positions.
  The model, i.e. the tagSNVs, is shown.
  Also annotations of tagSNVs
  like the minor allele is the base of another genome can
  be visualized.


}
\usage{

plotHaplotypeCluster(Lout,tagSNV,physPos=NULL,
   colRamp=12,val=c(0.0,2.0,1.0),chrom="",count=0,
   labelsNA=NULL,prange=NULL,labelsNA1 = c("model L") )


}
\arguments{
  \item{Lout}{the alleles of the individuals are provided; typically via a previous call
  of the \pkg{fabia} function \code{readSamplesSpfabia}.}
  \item{tagSNV}{tagSNVs given as numbers if all SNVs are enumerated.}
  \item{physPos}{physical position of tagSNVs in bp.}
  \item{colRamp}{passed to \code{snpPlot}: color representation.}
  \item{val}{vector of 3 components containing values for representation
    of reference allele, minor allele that is not tagSNV,
   minor allele that is tagSNV, respectively.}
  \item{chrom}{chromosome as string.}
  \item{count}{counter which is shown in the title if larger than 0 (for viewing a series of haplotype clusters).}
  \item{labelsNA}{labels for the individuals.}
  \item{prange}{vector of two integer values giving the begin and end of
   a subsegment of the segment for zooming in.}
  \item{labelsNA1}{labels for modelobtained by \pkg{fabia}; other tagSNV
    annotations like the agreement of the minor allele with other
    genomes can be provided.}


}
\details{

  A haplotype cluster is plotted showing tagSNVs and
  minor alleles of other SNVs.
  Provided are individuals to plot together with tagSNVs and their
  physical positions.
  Other annotations of tagSNVs can
  be visualized.

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

library(fabia)

data(hapRes)

data(simu)

namesL <- simu[["namesL"]]
haploN <- simu[["haploN"]]
snvs <- simu[["snvs"]]
annot <- simu[["annot"]]
alleleIimp <- simu[["alleleIimp"]]
write.table(namesL,file="dataSim1fabia_individuals.txt",
   quote = FALSE,row.names = FALSE,col.names = FALSE)
write(as.integer(haploN),file="dataSim1fabia_annot.txt",
   ncolumns=100)
write(as.integer(snvs),file="dataSim1fabia_annot.txt",
   append=TRUE,ncolumns=100)
write.table(annot,file="dataSim1fabia_annot.txt",
   sep = " ", quote = FALSE,row.names = FALSE,
   col.names = FALSE,append=TRUE)
write(as.integer(haploN),file="dataSim1fabia_mat.txt",
   ncolumns=100)
write(as.integer(snvs),file="dataSim1fabia_mat.txt",
   append=TRUE,ncolumns=100)

for (i in 1:haploN) {

  a1 <- which(alleleIimp[i,]>0.01)

  al <- length(a1)
  b1 <- alleleIimp[i,a1]

  a1 <- a1 - 1
  dim(a1) <- c(1,al)
  b1 <- format(as.double(b1),nsmall=1)
  dim(b1) <- c(1,al)

  write.table(al,file="dataSim1fabia_mat.txt",
     sep = " ", quote = FALSE,row.names = FALSE,
     col.names = FALSE,append=TRUE)
  write.table(a1,file="dataSim1fabia_mat.txt",
     sep = " ", quote = FALSE,row.names = FALSE,
     col.names = FALSE,append=TRUE)
  write.table(b1,file="dataSim1fabia_mat.txt",
     sep = " ", quote = FALSE,row.names = FALSE,
     col.names = FALSE,append=TRUE)

}



mergedHaploClusterList <- hapRes$mergedHaploClusterList
individuals <- individuals(mergedHaploClusterList[[1]])
tagSNVs <- tagSNVs(mergedHaploClusterList[[1]])
tagSNVs <-
   as.integer(sort.int(as.integer(unique(tagSNVs))))
tagSNVPositions <-
   tagSNVPositions(mergedHaploClusterList[[1]])
labelIndividuals <-
   labelIndividuals(mergedHaploClusterList[[1]])
Lout <- readSamplesSpfabia(X="dataSim1fabia_mat",
   samples=individuals,lowerB=0,upperB=1000.0)
tagSNVsL <- list(tagSNVs)
labelsK <- c("model L")
plotHaplotypeCluster(Lout=Lout,tagSNV=tagSNVsL,
   physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),
   chrom="1",count=0,labelsNA=labelIndividuals,
   labelsNA1=labelsK)

setwd(old_dir)

}



\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{hplot}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
