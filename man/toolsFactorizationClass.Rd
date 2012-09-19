\name{toolsFactorizationClass}

\alias{toolsFactorizationClass}
\alias{plotL}
\alias{plotL,Factorization-method}
\alias{plotL,Factorization,numeric,ANY,ANY,character,numeric,numeric,character,numeric-method}
\alias{histL}
\alias{histL,Factorization-method}
\alias{histL,Factorization,numeric,ANY,ANY,numeric,numeric-method}
\alias{topLZ}
\alias{topLZ,Factorization-method}
\alias{topLZ,Factorization,numeric,character,logical-method}
\alias{topLZ,Factorization,numeric,character,logical,ANY-method}
\alias{topLZ,Factorization,numeric,character,logical,ANY,ANY-method}



\title{Tools to analyze results of \pkg{fabia}}
\description{

  These tools allow to analyze results of the package \pkg{fabia}. They
  can be used to identify haplotype cluster regions and for adjusting
  the parameters of \code{extractHaploClusters} and \code{hapFabia} such as
  \code{ps} (top L values for extraction), \code{psZ} (top Z values for
  extraction), \code{inteA} (length of
  histogram bins).

  \code{plotL} plots the loadings of a \pkg{fabia} result that are above a threshold
  either as points, histogram or by a smooth scatter plot.

  \code{topLZ} returns largest L or Z values of a \pkg{fabia} result, where
  thresholds are given either by a quantile or by a value.

  \code{histL} supplies a histogram of the loadings
  obtained by \pkg{fabia}.

  


}
\usage{

plotL(res,n=1,p=NULL,w=NULL,type="points",
   intervv=500,off=0,t="p",cex=1)

histL(res,n=1,p=NULL,w=NULL,intervv=500,off=0)

topLZ(res,n=1,LZ="L",indices=TRUE,p=NULL,w=NULL)


}
\arguments{
  \item{res}{\pkg{fabia} result; instance of the class \code{Factorization}.}
  \item{n}{the number of the bicluster to consider.}
  \item{p}{the quantile threshold above which values are returned (p or w must be given).}
  \item{w}{the value threshold above which values are returned (p or w must be given).}
  \item{type}{the type of the plot: \code{type=c("points","histogram","smoothScatter")}.}
  \item{intervv}{length of the interval bins for histograms.}
  \item{off}{offset of the interval bins from zero for histograms.}
  \item{t}{points type for the plot.}
  \item{cex}{size of the points for the plot.}
  \item{LZ}{"L" for loadings L or "Z" for factors Z.}
  \item{indices}{if \code{TRUE} (default) indices are given and otherwise values.}

}
\details{


  \code{plotL} plots the loadings of a \pkg{fabia} result that are above a threshold
  either as points, histogram or by a smooth scatter plot.
  Thresholds can be given by a quantile or by a value.

  \code{topLZ} returns largest L or Z indices/values of a \pkg{fabia} result.
  Thresholds are given by quantile or by a value.

  \code{histL} computes histogram of the loadings obtained by \pkg{fabia}.

  Implementation in \R.

}

\value{
\item{}{\code{plotL}: nothing.}
\item{}{\code{topLZ}: vector of indices or values depending on the logical parameter indices.}
\item{}{\code{histL}: object of class \code{histogram}.}
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

data(res)

plotL(res,n=1,p=0.95,w=NULL,type="histogram",
   intervv=50,off=0,t="p",cex=1)
plotL(res,n=1,p=0.95,w=NULL,type="points",
   intervv=50,off=0,t="p",cex=1)
plotL(res,n=1,p=NULL,w=0.5,type="points",
   intervv=50,off=0,t="p",cex=1)
plotL(res,n=1,p=0.95,w=NULL,type="smooth",
   intervv=50,off=0,t="p",cex=1)
plotL(res,n=1,p=NULL,w=0.5,type="smooth",
   intervv=50,off=0,t="p",cex=1)

topLZ(res,n=1,LZ="L",indices=TRUE,p=0.95,w=NULL)
topLZ(res,n=1,LZ="L",indices=TRUE,p=NULL,w=0.95)

topLZ(res,n=1,LZ="Z",indices=TRUE,p=0.95,w=NULL)
topLZ(res,n=1,LZ="Z",indices=TRUE,p=NULL,w=0.4)

topLZ(res,n=1,LZ="L",indices=FALSE,p=0.95,w=NULL)
topLZ(res,n=1,LZ="L",indices=FALSE,p=NULL,w=0.95)

topLZ(res,n=1,LZ="Z",indices=FALSE,p=0.95,w=NULL)
topLZ(res,n=1,LZ="Z",indices=FALSE,p=NULL,w=0.4)

h1 <- histL(res,n=1,p=0.9,w=NULL,intervv=50,off=0)
print(h1$counts)
h1 <- histL(res,n=1,p=NULL,w=0.5,intervv=50,off=0)
print(h1$counts)

}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
