\name{plotIBDsegment}

\alias{plotIBDsegment}

\title{Plots an IBD segment given genotype data and tagSNVs}

\description{

  \code{plotIBDsegment}: \R implementation of \code{plotIBDsegment}.

  A IBD segment is plotted where individuals that
  are plotted must be provided together with tagSNVs and their
  physical positions. The individuals are provided as genotyping matrix.
  The model, i.e. the tagSNVs, is shown, too.
  Annotations of tagSNVs
  like match with another genome (Neandertal, Denisova) can
  be visualized.


}
\usage{

plotIBDsegment(Lout,tagSNV,physPos=NULL,colRamp=12,val=c(0.0,2.0,1.0),chrom="",count=0,labelsNA=NULL,prange=NULL,labelsNA1=c("model L"),grid=FALSE,pairs=FALSE,...)


}
\arguments{
  \item{Lout}{the alleles of the individuals are provided; typically via a previous call
  of the \pkg{fabia} function \code{readSamplesSpfabia}.}
  \item{tagSNV}{tagSNVs given as numbers if all SNVs are enumerated.}
  \item{physPos}{physical position of tagSNVs in bp.}
  \item{colRamp}{passed to \code{matrixPlot}: color representation.}
  \item{val}{vector of 3 components containing values for representation
    of reference allele, minor allele that is not tagSNV,
   minor allele that is tagSNV, respectively.}
  \item{chrom}{chromosome as string.}
  \item{count}{counter which is shown in the title if larger than 0 (for viewing a series of IBD segments).}
  \item{labelsNA}{labels for the individuals.}
  \item{prange}{vector of two integer values giving the begin and end of
   a subinterval of the interval for zooming in.}
 \item{labelsNA1}{labels for tagSNVs by the model obtained by \pkg{fabia}; other
   tagSNV annotations like match with another genome (Neandertal, Denisova) can
  be visualized.}
  \item{grid}{does the plot have a grid?; default FALSE (no).}
  \item{pairs}{for pairwise groups, e.g. case-control, twins, etc.; default FALSE (no).}
  \item{...}{other graphical parameters may also be passed as arguments to this function.}


}
\details{

  A IBD segment is plotted showing tagSNVs and
  minor alleles of other SNVs.
  Provided are individuals to plot together with tagSNVs and their
  physical positions.
  Other annotations of tagSNVs can
  be visualized.

In the plot the $y$-axis
gives the individuals or the chromosomes
and the $x$-axis consecutive SNVs. The default color coding uses yellow
for major alleles, violet for minor alleles of tagSNVs, and blue
for minor alleles of other SNVs.
\code{model L} indicates tagSNVs
identified by hapFabia in violet.

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



mergedIBDsegmentList <- hapRes$mergedIBDsegmentList
individuals <- individuals(mergedIBDsegmentList[[1]])
tagSNVs <- tagSNVs(mergedIBDsegmentList[[1]])
tagSNVs <-
   as.integer(sort.int(as.integer(unique(tagSNVs))))
tagSNVPositions <-
   tagSNVPositions(mergedIBDsegmentList[[1]])
labelIndividuals <-
   labelIndividuals(mergedIBDsegmentList[[1]])
Lout <- readSamplesSpfabia(X="dataSim1fabia_mat",
   samples=individuals,lowerB=0,upperB=1000.0)
tagSNVsL <- list(tagSNVs)
labelsK <- c("model L")
plotIBDsegment(Lout=Lout,tagSNV=tagSNVsL,
   physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),
   chrom="1",count=0,labelsNA=labelIndividuals,
   labelsNA1=labelsK)

setwd(old_dir)

}
\value{Plots an IBD segment given genotype data and tagSNVs}


\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{hplot}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
