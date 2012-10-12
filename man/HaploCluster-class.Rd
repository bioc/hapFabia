\name{HaploCluster-class}
\docType{class}

% Class
\alias{HaploCluster-class}
\alias{HaploCluster-method}

% Constructor-like function:
\alias{HaploCluster}
\alias{HaploCluster,ANY-method}
\alias{HaploCluster,HaploCluster-method}
\alias{HaploCluster,numeric,numeric,character,numeric,numeric,numeric,numeric,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector,vector-method}


% Coercion methods:

% Accessor methods:
\alias{ID}
\alias{ID,HaploCluster-method}
\alias{ID<-}
\alias{ID<-,HaploCluster,numeric-method}

\alias{bicluster_id}
\alias{bicluster_id,HaploCluster-method}
\alias{bicluster_id<-}
\alias{bicluster_id<-,HaploCluster,numeric-method}

\alias{chromosome}
\alias{chromosome,HaploCluster-method}
\alias{chromosome<-}
\alias{chromosome<-,HaploCluster,character-method}

\alias{haploClusterPos}
\alias{haploClusterPos,HaploCluster-method}
\alias{haploClusterPos<-}
\alias{haploClusterPos<-,HaploCluster,numeric-method}

\alias{haploClusterLength}
\alias{haploClusterLength,HaploCluster-method}
\alias{haploClusterLength<-}
\alias{haploClusterLength<-,HaploCluster,numeric-method}

\alias{numberIndividuals}
\alias{numberIndividuals,HaploCluster-method}
\alias{numberIndividuals<-}
\alias{numberIndividuals<-,HaploCluster,numeric-method}

\alias{numbertagSNVs}
\alias{numbertagSNVs,HaploCluster-method}
\alias{numbertagSNVs<-}
\alias{numbertagSNVs<-,HaploCluster,numeric-method}

\alias{individuals}
\alias{individuals,HaploCluster-method}
\alias{individuals<-}
\alias{individuals<-,HaploCluster,vector-method}

\alias{tagSNVs}
\alias{tagSNVs,HaploCluster-method}
\alias{tagSNVs<-}
\alias{tagSNVs<-,HaploCluster,vector-method}

\alias{populationIndividuals}
\alias{populationIndividuals,HaploCluster-method}
\alias{populationIndividuals<-}
\alias{populationIndividuals<-,HaploCluster,vector-method}

\alias{idIndividuals}
\alias{idIndividuals,HaploCluster-method}
\alias{idIndividuals<-}
\alias{idIndividuals<-,HaploCluster,vector-method}

\alias{labelIndividuals}
\alias{labelIndividuals,HaploCluster-method}
\alias{labelIndividuals<-}
\alias{labelIndividuals<-,HaploCluster,vector-method}

\alias{platformIndividuals}
\alias{platformIndividuals,HaploCluster-method}
\alias{platformIndividuals<-}
\alias{platformIndividuals<-,HaploCluster,vector-method}

\alias{coreClusterIndividuals}
\alias{coreClusterIndividuals,HaploCluster-method}
\alias{coreClusterIndividuals<-}
\alias{coreClusterIndividuals<-,HaploCluster,vector-method}

\alias{tagSNVPositions}
\alias{tagSNVPositions,HaploCluster-method}
\alias{tagSNVPositions<-}
\alias{tagSNVPositions<-,HaploCluster,vector-method}

\alias{tagSNVAlleles}
\alias{tagSNVAlleles,HaploCluster-method}
\alias{tagSNVAlleles<-}
\alias{tagSNVAlleles<-,HaploCluster,vector-method}

\alias{tagSNVNames}
\alias{tagSNVNames,HaploCluster-method}
\alias{tagSNVNames<-}
\alias{tagSNVNames<-,HaploCluster,vector-method}

\alias{tagSNVFreq}
\alias{tagSNVFreq,HaploCluster-method}
\alias{tagSNVFreq<-}
\alias{tagSNVFreq<-,HaploCluster,vector-method}

\alias{tagSNVGroupFreq}
\alias{tagSNVGroupFreq,HaploCluster-method}
\alias{tagSNVGroupFreq<-}
\alias{tagSNVGroupFreq<-,HaploCluster,vector-method}

\alias{tagSNVChange}
\alias{tagSNVChange,HaploCluster-method}
\alias{tagSNVChange<-}
\alias{tagSNVChange<-,HaploCluster,vector-method}

\alias{tagSNVsPerIndividual}
\alias{tagSNVsPerIndividual,HaploCluster-method}
\alias{tagSNVsPerIndividual<-}
\alias{tagSNVsPerIndividual<-,HaploCluster,vector-method}

\alias{individualPerTagSNV}
\alias{individualPerTagSNV,HaploCluster-method}
\alias{individualPerTagSNV<-}
\alias{individualPerTagSNV<-,HaploCluster,vector-method}

\alias{tagSNVAnno}
\alias{tagSNVAnno,HaploCluster-method}
\alias{tagSNVAnno<-}
\alias{tagSNVAnno<-,HaploCluster,vector-method}



% summary method:
\alias{summary,HaploCluster-method}
% plot method:
\alias{plot,HaploCluster-method}
\alias{plot,HaploCluster,missing-method}


\title{HaploCluster instances and methods}

\description{\code{HaploCluster} is a class to store characteristics of
  a haplotype cluster in one of its instances.
  Characteristics of a haplotype cluster include its genomic
  position, its length, the individuals/chromosomes that belong to it,
  the tagSNVs that tag/mark it, etc.
}

\section{Slots}{
  Objects of class \code{HaploCluster} have the following slots:
  \describe{
    
  \item{\code{ID}}{number of the haplotype cluster in the current extraction.}
  \item{\code{bicluster_id}}{ID of the bicluster the haplotype cluster was found in.}
  \item{\code{chromosome}}{the chromosome.}
  \item{\code{haploClusterPos}}{genomic location of the haplotype cluster.}
  \item{\code{haploClusterLength}}{length of the haplotype cluster in
  the number of SNVs. For the length in kbp: \code{max(tagSNVPositions(x))-min(tagSNVPositions(x))}.}
  \item{\code{numberIndividuals}}{number of samples belonging to the haplotype cluster.}
  \item{\code{numbertagSNVs}}{number tagSNVs marking the haplotype cluster.}
  \item{\code{individuals}}{IDs of individuals or chromosomes belonging to the haplotype cluster.}
  \item{\code{tagSNVs}}{IDs of SNVs that mark the haplotype cluster (tagSNVs).}
  \item{\code{populationIndividuals}}{the population each individual belongs to.}
  \item{\code{idIndividuals}}{IDs of the individuals or chromosomes.}
  \item{\code{labelIndividuals}}{label of the individuals.}
  \item{\code{platformIndividuals}}{for each individual the technology/platform that was used to genotype it.}
  \item{\code{coreClusterIndividuals}}{IDs of individuals that constitute the core
  of the haplotype cluster.}
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
 Constructor of class HaploCluster.
 
  \code{HaploCluster(ID=0,bicluster_id=0,chromosome="",haploClusterPos=0,haploClusterLength=0,numberIndividuals=0,numbertagSNVs=0,individuals=as.vector(0),tagSNVs=as.vector(0),populationIndividuals=as.vector(""),idIndividuals=as.vector(0),labelIndividuals=as.vector(""),platformIndividuals=as.vector(""),coreClusterIndividuals=as.vector(0),tagSNVPositions=as.vector(0),tagSNVAlleles=as.vector(""),tagSNVNames=as.vector(""),tagSNVFreq=as.vector(0),tagSNVGroupFreq=as.vector(0),tagSNVChange=as.vector(0),tagSNVsPerIndividual=as.vector(0),individualPerTagSNV=as.vector(0),tagSNVAnno=as.vector(c(list(""))))}
}

\section{Accessors}{
  In the following \code{x} denotes a HaploCluster object.
   
 
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
      \code{haploClusterPos(x)}, \code{haploClusterPos(x) <- value}:
      Returns or sets \code{haploClusterPos}, where the return value and
      \code{value} are both numeric.
    }



    \item{}{
      \code{haploClusterLength(x)}, \code{haploClusterLength(x) <- value}:
      Returns or sets \code{haploClusterLength}, where the return value and
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


\S4method{plot}{HaploCluster}(x,filename, ...)

\S4method{summary}{HaploCluster}(object, ...) 

}


\section{Signatures}{
  \describe{


\item{plot}{\code{signature(x = "HaploCluster", y = "missing")}}{ Plot
  of a haplotype cluster, where tagSNVs, minor and major alleles are plotted.}


\item{summary}{\code{signature(object = "HaploCluster")}}{ Summary of
  haplotype cluster object.}


}
}

\arguments{
  \item{x}{object of the class \code{HaploCluster}.}
  \item{object}{object of the class \code{HaploCluster}.}
  \item{filename}{filename of the file that contains the genotyping
    data. Call of \code{readSamplesSpfabia} from the package \pkg{fabia}
  to read genotype data for the plot. ATTENTION: without file type ".txt"!}
  \item{...}{further arguments.}

}


\value{
  \item{}{no value.}
}




\details{

  \describe{

    \item{\bold{Plot}}{Plots a haplotype cluster where the SNVs within
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



\item{\bold{Summary}}{Prints the ID of the haplotype cluster, the
  bicluster it is derived from, the chromosome and position where it is
  located, its length, the number of individuals/chromosomes that belong
  to it, the number of tagSNVs that mark it.}

}  

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

haploCluster <- mergedHaploClusterList[[1]]


# Summary method
summary(haploCluster)


# Plot method
plot(haploCluster,filename="dataSim1fabia_mat")



# ACCESSORS 

# IDs of the haplotype cluster
   ID(haploCluster)
   bicluster_id(haploCluster)

# General Information
    haploClusterPos(haploCluster)
    haploClusterLength(haploCluster)
    numberIndividuals(haploCluster)
    numbertagSNVs(haploCluster)
    coreClusterIndividuals(haploCluster)

# Information on individuals / chromosomes
    individuals(haploCluster)
    populationIndividuals(haploCluster)
    idIndividuals(haploCluster)
    labelIndividuals(haploCluster)
    platformIndividuals(haploCluster)
    tagSNVsPerIndividual(haploCluster)

# Information on tagSNVs
    tagSNVs(haploCluster)
    tagSNVPositions(haploCluster)
    tagSNVAlleles(haploCluster)
    tagSNVNames(haploCluster)
    tagSNVFreq(haploCluster)
    tagSNVGroupFreq(haploCluster)
    tagSNVChange(haploCluster)
    individualPerTagSNV(haploCluster)
    tagSNVAnno(haploCluster)


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
\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
