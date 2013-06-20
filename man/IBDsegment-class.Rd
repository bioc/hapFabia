\name{IBDsegment-class}
\docType{class}

% Class
\alias{IBDsegment-class}
\alias{IBDsegment-method}

% Constructor-like function:
\alias{IBDsegment}
\alias{IBDsegment,ANY-method}
\alias{IBDsegment,IBDsegment-method}
\alias{IBDsegment,numeric,numeric,character,numeric,numeric,numeric,numeric,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector-method}


% Coercion methods:

% Accessor methods:
\alias{ID}
\alias{ID,IBDsegment-method}
\alias{ID<-}
\alias{ID<-,IBDsegment,numeric-method}

\alias{bicluster_id}
\alias{bicluster_id,IBDsegment-method}
\alias{bicluster_id<-}
\alias{bicluster_id<-,IBDsegment,numeric-method}

\alias{chromosome}
\alias{chromosome,IBDsegment-method}
\alias{chromosome<-}
\alias{chromosome<-,IBDsegment,character-method}

\alias{IBDsegmentPos}
\alias{IBDsegmentPos,IBDsegment-method}
\alias{IBDsegmentPos<-}
\alias{IBDsegmentPos<-,IBDsegment,numeric-method}

\alias{IBDsegmentLength}
\alias{IBDsegmentLength,IBDsegment-method}
\alias{IBDsegmentLength<-}
\alias{IBDsegmentLength<-,IBDsegment,numeric-method}

\alias{numberIndividuals}
\alias{numberIndividuals,IBDsegment-method}
\alias{numberIndividuals<-}
\alias{numberIndividuals<-,IBDsegment,numeric-method}

\alias{numbertagSNVs}
\alias{numbertagSNVs,IBDsegment-method}
\alias{numbertagSNVs<-}
\alias{numbertagSNVs<-,IBDsegment,numeric-method}

\alias{individuals}
\alias{individuals,IBDsegment-method}
\alias{individuals<-}
\alias{individuals<-,IBDsegment,vector-method}

\alias{tagSNVs}
\alias{tagSNVs,IBDsegment-method}
\alias{tagSNVs<-}
\alias{tagSNVs<-,IBDsegment,vector-method}

\alias{populationIndividuals}
\alias{populationIndividuals,IBDsegment-method}
\alias{populationIndividuals<-}
\alias{populationIndividuals<-,IBDsegment,vector-method}

\alias{idIndividuals}
\alias{idIndividuals,IBDsegment-method}
\alias{idIndividuals<-}
\alias{idIndividuals<-,IBDsegment,vector-method}

\alias{labelIndividuals}
\alias{labelIndividuals,IBDsegment-method}
\alias{labelIndividuals<-}
\alias{labelIndividuals<-,IBDsegment,vector-method}

\alias{platformIndividuals}
\alias{platformIndividuals,IBDsegment-method}
\alias{platformIndividuals<-}
\alias{platformIndividuals<-,IBDsegment,vector-method}

\alias{coreClusterIndividuals}
\alias{coreClusterIndividuals,IBDsegment-method}
\alias{coreClusterIndividuals<-}
\alias{coreClusterIndividuals<-,IBDsegment,vector-method}

\alias{tagSNVPositions}
\alias{tagSNVPositions,IBDsegment-method}
\alias{tagSNVPositions<-}
\alias{tagSNVPositions<-,IBDsegment,vector-method}

\alias{tagSNVAlleles}
\alias{tagSNVAlleles,IBDsegment-method}
\alias{tagSNVAlleles<-}
\alias{tagSNVAlleles<-,IBDsegment,vector-method}

\alias{tagSNVNames}
\alias{tagSNVNames,IBDsegment-method}
\alias{tagSNVNames<-}
\alias{tagSNVNames<-,IBDsegment,vector-method}

\alias{tagSNVFreq}
\alias{tagSNVFreq,IBDsegment-method}
\alias{tagSNVFreq<-}
\alias{tagSNVFreq<-,IBDsegment,vector-method}

\alias{tagSNVGroupFreq}
\alias{tagSNVGroupFreq,IBDsegment-method}
\alias{tagSNVGroupFreq<-}
\alias{tagSNVGroupFreq<-,IBDsegment,vector-method}

\alias{tagSNVChange}
\alias{tagSNVChange,IBDsegment-method}
\alias{tagSNVChange<-}
\alias{tagSNVChange<-,IBDsegment,vector-method}

\alias{tagSNVsPerIndividual}
\alias{tagSNVsPerIndividual,IBDsegment-method}
\alias{tagSNVsPerIndividual<-}
\alias{tagSNVsPerIndividual<-,IBDsegment,vector-method}

\alias{individualPerTagSNV}
\alias{individualPerTagSNV,IBDsegment-method}
\alias{individualPerTagSNV<-}
\alias{individualPerTagSNV<-,IBDsegment,vector-method}

\alias{tagSNVAnno}
\alias{tagSNVAnno,IBDsegment-method}
\alias{tagSNVAnno<-}
\alias{tagSNVAnno<-,IBDsegment,vector-method}



% summary method:
\alias{summary,IBDsegment-method}
% plot method:
\alias{plot,IBDsegment-method}
\alias{plot,IBDsegment,missing-method}
% extended plot method:
\alias{plotLarger}
\alias{plotLarger,IBDsegment-method}
\alias{plotLarger,IBDsegment,character,numeric-method}
\alias{plotLarger,IBDsegment,missing-method}


\title{IBDsegment instances and methods}

\description{\code{IBDsegment} is a class to store characteristics of
  an IBD segment in one of its instances.
  Characteristics of an IBD segment include its genomic
  position, its length, the individuals/chromosomes that belong to it,
  the tagSNVs that tag/mark it, etc.
}

\section{Slots}{
  Objects of class \code{IBDsegment} have the following slots:
  \describe{
    
  \item{\code{ID}}{number of the IBD segment in the current extraction.}
  \item{\code{bicluster_id}}{ID of the bicluster the IBD segment was found in.}
  \item{\code{chromosome}}{the chromosome.}
  \item{\code{IBDsegmentPos}}{genomic location of the IBD segment.}
  \item{\code{IBDsegmentLength}}{length of the IBD segment in
  the number of SNVs. For the length in bp: \code{max(tagSNVPositions(x))-min(tagSNVPositions(x))}.}
  \item{\code{numberIndividuals}}{number of samples belonging to the IBD segment.}
  \item{\code{numbertagSNVs}}{number tagSNVs marking the IBD segment.}
  \item{\code{individuals}}{IDs of individuals or chromosomes belonging to the IBD segment.}
  \item{\code{tagSNVs}}{IDs of SNVs that mark the IBD segment (tagSNVs).}
  \item{\code{populationIndividuals}}{the population each individual belongs to.}
  \item{\code{idIndividuals}}{IDs of the individuals or chromosomes.}
  \item{\code{labelIndividuals}}{label of the individuals.}
  \item{\code{platformIndividuals}}{for each individual the technology/platform that was used to genotype it.}
  \item{\code{coreClusterIndividuals}}{IDs of individuals that constitute the core
  of the IBD segment.}
  \item{\code{tagSNVPositions}}{physical positions of the tagSNVs on the
  chromosome in base pairs.}
  \item{\code{tagSNVAlleles}}{alleles of the tagSNVs in the form \code{Ref:Alt} where Ref denotes reference allele and Alt the alternative allele.}
  \item{\code{tagSNVNames}}{name of the tagSNVs according to a given annotation.}
  \item{\code{tagSNVFreq}}{frequency of tagSNVs in the whole data set.}
  \item{\code{tagSNVGroupFreq}}{frequency of tagSNVs in the population that is considered.}
  \item{\code{tagSNVChange}}{if the minor allele was more frequent than the
  major, then both were switched. Switching is marked by a 1 and
  otherwise it is 0.}
  \item{\code{tagSNVsPerIndividual}}{for each sample: tagSNVs are counted for which the sample has the minor allele.}
  \item{\code{individualPerTagSNV}}{for each tagSNV: samples are counted for which the SNV has its minor allele.}
\item{\code{tagSNVAnno}}{the functional annotation of tagSNVs for each tagSNV: like stop-loss,
  stop-gain, non-synonymous, synonymous, promoter, exonic, intronic, intergenic, etc.}

     }
  }

\section{Constructor}{
 Constructor of class IBDsegment.
 
  \code{IBDsegment(ID=0,bicluster_id=0,chromosome="",IBDsegmentPos=0,IBDsegmentLength=0,numberIndividuals=0,numbertagSNVs=0,individuals=as.vector(0),tagSNVs=as.vector(0),populationIndividuals=as.vector(""),idIndividuals=as.vector(0),labelIndividuals=as.vector(""),platformIndividuals=as.vector(""),coreClusterIndividuals=as.vector(0),tagSNVPositions=as.vector(0),tagSNVAlleles=as.vector(""),tagSNVNames=as.vector(""),tagSNVFreq=as.vector(0),tagSNVGroupFreq=as.vector(0),tagSNVChange=as.vector(0),tagSNVsPerIndividual=as.vector(0),individualPerTagSNV=as.vector(0),tagSNVAnno=as.vector(c(list(""))))}
}

\section{Accessors}{
  In the following \code{x} denotes an IBDsegment object.
   
 
  \describe{

    \item{}{
      \code{ID(x)}, \code{ID(x) <- value}:
      Returns or sets \code{ID}, where the return value and
      \code{value} are both numeric.
    }

    \item{}{
      \code{bicluster_id(x)}, \code{bicluster_id(x) <- value}:
      Returns or sets \code{bicluster_id}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{chromosome(x)}, \code{chromosome(x) <- value}:
      Returns or sets \code{chromosome}, where the return value and
      \code{value} are both strings.
    }



    \item{}{
      \code{IBDsegmentPos(x)}, \code{IBDsegmentPos(x) <- value}:
      Returns or sets \code{IBDsegmentPos}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{IBDsegmentLength(x)}, \code{IBDsegmentLength(x) <- value}:
      Returns or sets \code{IBDsegmentLength}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{numberIndividuals(x)}, \code{numberIndividuals(x) <- value}:
      Returns or sets \code{numberIndividuals}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{numbertagSNVs(x)}, \code{numbertagSNVs(x) <- value}:
      Returns or sets \code{numbertagSNVs}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{individuals(x)}, \code{individuals(x) <- value}:
      Returns or sets \code{individuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVs(x)}, \code{tagSNVs(x) <- value}:
      Returns or sets \code{tagSNVs}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{populationIndividuals(x)}, \code{populationIndividuals(x) <- value}:
      Returns or sets \code{populationIndividuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{idIndividuals(x)}, \code{idIndividuals(x) <- value}:
      Returns or sets \code{idIndividuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{labelIndividuals(x)}, \code{labelIndividuals(x) <- value}:
      Returns or sets \code{labelIndividuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{platformIndividuals(x)}, \code{platformIndividuals(x) <- value}:
      Returns or sets \code{platformIndividuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{coreClusterIndividuals(x)}, \code{coreClusterIndividuals(x) <- value}:
      Returns or sets \code{coreClusterIndividuals}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVPositions(x)}, \code{tagSNVPositions(x) <- value}:
      Returns or sets \code{tagSNVPositions}, where the return value and
      \code{value} are both vectors. 
    }



    \item{}{
      \code{tagSNVAlleles(x)}, \code{tagSNVAlleles(x) <- value}:
      Returns or sets \code{tagSNVAlleles}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVNames(x)}, \code{tagSNVNames(x) <- value}:
      Returns or sets \code{tagSNVNames}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVFreq(x)}, \code{tagSNVFreq(x) <- value}:
      Returns or sets \code{tagSNVFreq}, where the return value and
      \code{value} are both vectors. 
    }



    \item{}{
      \code{tagSNVGroupFreq(x)}, \code{tagSNVGroupFreq(x) <- value}:
      Returns or sets \code{tagSNVGroupFreq}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVChange(x)}, \code{tagSNVChange(x) <- value}:
      Returns or sets \code{tagSNVChange}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVsPerIndividual(x)}, \code{tagSNVsPerIndividual(x) <- value}:
      Returns or sets \code{tagSNVsPerIndividual}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{individualPerTagSNV(x)}, \code{individualPerTagSNV(x) <- value}:
      Returns or sets \code{individualPerTagSNV}, where the return value and
      \code{value} are both vectors.
    }



    \item{}{
      \code{tagSNVAnno(x)}, \code{tagSNVAnno(x) <- value}:
      Returns or sets \code{tagSNVAnno}, where the return value and
      \code{value} are both vectors.
    }




  }
}


\usage{


\S4method{plot}{IBDsegment}(x,filename, ...)

\S4method{plotLarger}{IBDsegment}(x,filename,fact=1.0,addSamp=c(), ...)

\S4method{summary}{IBDsegment}(object, ...) 

}


\section{Signatures}{
  \describe{


\item{plot}{\code{signature(x = "IBDsegment", y = "missing")}}{ Plot
  of an IBD segment, where tagSNVs, minor and major alleles are plotted.}

\item{plotLarger}{\code{signature(x="IBDsegment", filename="character",fact="numeric",addSamp="ANY")}}{ Plot
  of an IBD segment with additional individuals and the IBD segment
  extended to the left and to the right.}


\item{summary}{\code{signature(object = "IBDsegment")}}{ Summary of
  IBD segment object.}


}
}

\arguments{
  \item{x}{object of the class \code{IBDsegment}.}
  \item{object}{object of the class \code{IBDsegment}.}
  \item{filename}{filename of the file that contains the genotyping
    data. Call of \code{readSamplesSpfabia} from the package \pkg{fabia}
  to read genotype data for the plot. ATTENTION: without file type ".txt"!}
  \item{fact}{factor by which the IBD segment is extended. The extension
  is done both at the left and at the right hand side.}
  \item{addSamp}{vector giving additional individuals, e.g. for adding
  10 individuals out of N the code \code{sample(N,10)} may be used.}
  \item{...}{further arguments.}

}


\value{
  \item{}{no value.}
}




\details{

  \describe{

    \item{\bold{plot}}{Plots an IBD segment where the SNVs within
    the cluster are shown. In the plot the $y$-axis
gives the individuals or the chromosomes
and the $x$-axis consecutive SNVs.
Minor alleles of tagSNVs are marked by a
    particular color. Minor
    alleles of other SNVs (private, common, from other clusters, etc) are
    marked by a different color.  The model from \pkg{fabia} is also
    shown.
     The default color coding uses yellow
for major alleles, violet for minor alleles of tagSNVs, and blue
for minor alleles of other SNVs.
\code{model L} indicates tagSNVs
identified by hapFabia in violet.
}

    \item{\bold{plotLarger}}{Plots an IBD segment where the SNVs within
    the cluster are shown, however additional individuals that do not
    possess the IBD segment can be added.
    Further the range can be increased by a certain factor.
    Additional arguments are \code{fact} and \code{addSamp}.
    \code{fact} gives the factor by which the IBD segment should be
    extended to both the left and the right.
    \code{addSamp} is a vector of individuals/haplotypes
    which are added to the individuals that possess the IBD segment.
    The plot allows to view the IBD segment in the context of its
    DNA location and additional individuals. 
    }



\item{\bold{summary}}{Prints the ID of the IBD segment, the
  bicluster it is derived from, the chromosome and position where it is
  located, its length, the number of individuals/chromosomes that belong
  to it, the number of tagSNVs that mark it.}

}  

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

IBDsegment <- mergedIBDsegmentList[[1]]


# Summary method
summary(IBDsegment)


# Plot method
plot(IBDsegment,filename="dataSim1fabia_mat")

# Extended plot: more examples and borders 
plotLarger(IBDsegment,filename="dataSim1fabia_mat",3,sample(100,10))


# ACCESSORS 

# IDs of the IBD segment
   ID(IBDsegment)
   bicluster_id(IBDsegment)

# General Information
    IBDsegmentPos(IBDsegment)
    IBDsegmentLength(IBDsegment)
    numberIndividuals(IBDsegment)
    numbertagSNVs(IBDsegment)
    coreClusterIndividuals(IBDsegment)

# Information on individuals / chromosomes
    individuals(IBDsegment)
    populationIndividuals(IBDsegment)
    idIndividuals(IBDsegment)
    labelIndividuals(IBDsegment)
    platformIndividuals(IBDsegment)
    tagSNVsPerIndividual(IBDsegment)

# Information on tagSNVs
    tagSNVs(IBDsegment)
    tagSNVPositions(IBDsegment)
    tagSNVAlleles(IBDsegment)
    tagSNVNames(IBDsegment)
    tagSNVFreq(IBDsegment)
    tagSNVGroupFreq(IBDsegment)
    tagSNVChange(IBDsegment)
    individualPerTagSNV(IBDsegment)
    tagSNVAnno(IBDsegment)


setwd(old_dir)


}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}

\keyword{classes}
\keyword{methods}
\keyword{hplot}
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
