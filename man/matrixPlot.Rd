\name{matrixPlot}

\alias{matrixPlot}

\title{Basic plot function for IBD segments}

\description{

  \code{matrixPlot}: \R implementation of \code{matrixPlot}.

  Plots a matrix where different values are coded by different colors.
  Basically the image plot function \code{image} with a
  particular scaling, color coding, and axis.

}
\usage{

matrixPlot(x,range=NULL,yLabels=NULL,zlim=NULL,title=NULL,colRamp=12,grid=FALSE,pairs=FALSE,padj=NA,...)


}
\arguments{
  \item{x}{matrix that codes alleles and annotations of an IBD segment.}
  \item{range}{optional: physical range of the IBD segment.}
  \item{yLabels}{optional: labels of the individuals.}
  \item{zlim}{optional: limits imposed onto the matrix values.}
  \item{title}{title of the plot.}
  \item{colRamp}{color representation.}
  \item{grid}{does the plot have a grid?; default FALSE (no).}
  \item{pairs}{for pairwise groups, e.g. case-control, twins, etc.; default FALSE (no).}
  \item{padj}{adjustment for each tick label perpendicular to the reading direction.}
  \item{...}{other graphical parameters may also be passed as arguments to this function.}
}
\details{

 

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

mat <- matrix(0,nrow=10,ncol=40)
v1 <- sample(1:10,5)
v21 <- sample(1:40,4)
v22 <- sample(1:40,4)
w1 <- rep(0,10)
w2 <- rep(0,40)
w1[v1] <- 1
w2[v21] <- 1
w2[v22] <- 2
mat <- mat + tcrossprod(w1,w2)

matrixPlot(mat)


}
\value{Plots a matrix where different values are coded by different colors}
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{hplot}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
