\name{HaploClusterList-class}
\docType{class}

% Class
\alias{HaploClusterList-class}
\alias{HaploClusterList-method}

% Constructor-like function:
\alias{HaploClusterList}
\alias{HaploClusterList,ANY-method}
\alias{HaploClusterList,list,numeric,list-method}

% Coercion methods:

% Accessor methods:

\alias{haploClusters}
\alias{haploClusters,HaploClusterList-method}
\alias{haploClusters<-}
\alias{haploClusters<-,HaploClusterList,list-method}


\alias{lengthList}
\alias{lengthList,HaploClusterList-method}
\alias{lengthList<-}
\alias{lengthList<-,HaploClusterList,numeric-method}


\alias{statistics}
\alias{statistics,HaploClusterList-method}
\alias{statistics<-}
\alias{statistics<-,HaploClusterList,list-method}


\alias{[[}
\alias{[[,HaploClusterList-method}
\alias{[[,HaploClusterList,numeric,missing-method}
\alias{[[<-}
\alias{[[<-,HaploClusterList-method}
\alias{[[<-,HaploClusterList,numeric,missing,HaploCluster-method}


\alias{[}
\alias{[,HaploClusterList-method}
\alias{[,HaploClusterList,numeric,missing-method}
\alias{[<-}
\alias{[<-,HaploClusterList-method}
\alias{[<-,HaploClusterList,numeric,missing,HaploClusterList-method}



% summary method:
\alias{summary,HaploClusterList-method}
% plot method:
\alias{plot,HaploClusterList-method}
\alias{plot,HaploClusterList,missing-method}



\title{HaploClusterList instances and methods}

\description{\code{HaploClusterList} is a class to store a list of
  haplotype clusters with its statistics.
  Lists can be merged or analyzed in subsequent steps.
}

\section{Slots}{
  Objects of class \code{HaploClusterList} have the following slots:
  \describe{
        \item{\code{haploClusters}:}{List of haplotype clusters.}
         \item{\code{lengthList}:}{Number of haplotype clusters in the list.}
         \item{\code{statistics}:}{Statistics of haplotype clusters like
  average length, average number of individuals belonging to a haplotype
  cluster, average number of tagSNVs of a haplotype cluster, etc.}

     }
  }

\section{Constructor}{
 Constructor of class HaploClusterList.
 
 \code{HaploClusterList(haploClusters=list(),lengthList=0,statistics=list())}
   
}

\section{Accessors}{
  In the following \code{x} denotes a HaploClusterList object.
   
 
  \describe{

    \item{}{
      \code{haploClusters(x)}, \code{haploClusters(x) <- value}:
      Returns or sets \code{haploClusters}, where the return value and
      \code{value} are both a list. 
    }

    \item{}{
      \code{lengthList(x)}, \code{lengthList(x) <- value}:
      Returns or sets \code{lengthList}, where the return value and
      \code{value} are both a number. 
    }

    \item{}{
      \code{statistics(x)}, \code{statistics(x) <- value}:
      Returns or sets \code{statistics}, where the return value and
      \code{value} are both a list. 
    }


    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Returns or sets an entry in the list \code{x},
      where the return value and
      \code{value} are both an instance of the class \code{HaploCluster}. 
    }

    \item{}{
      \code{x[i]}, \code{x[i] <- value}:
      Returns or sets a sublist of the list \code{x},
      where the return value and
      \code{value} are both an instance of the class \code{HaploClusterList}. 
    }


  }
}


\usage{


\S4method{plot}{HaploClusterList}(x, ... )


\S4method{summary}{HaploClusterList}(object, ...) 

}


\section{Signatures}{
  \describe{


\item{plot}{\code{signature(x = "HaploClusterList", y = "missing")}}{
  Plotting all haplotype clusters of the list using an interactive command.}


\item{summary}{\code{signature(object = "HaploClusterList")}}{ Summary
  of a list of haplotype clusters, where the number of clusters and a
  statistics are given. }


}
}

\arguments{
  \item{x}{object of the class \code{HaploCluster}.}
  \item{object}{object of the class \code{HaploCluster}.}
  \item{...}{further arguments. PLOT: \code{filename}
    is the name of the file that contains the genotyping
    data used for plotting.
    The name is required for the call of \code{readSamplesSpfabia}
    from the package \pkg{fabia}
    to read genotype data for the plot. ATTENTION: \code{filename}
    without file type ".txt"!}
}


\value{
  \item{}{no value.}
}



\section{Functions that return objects of this class}{
  HaploClusterList objects are returned by \code{fabia}, \code{fabias}, \code{fabiap},
  \code{fabiasp},  \code{mfsc}, \code{nmfsc},
  \code{nmfdiv}, and \code{nmfeu}.} 

\section{Extension to store results of other methods}{
  
  The class \code{HaploClusterList} may contain the result of different matrix factorization
  methods. The methods may be generative or not.

  Methods my be "singular value decomposition" (M contains singular values
  as well as avini, L and Z are orthonormal matrices),
  "independent component analysis" (Z contains the projection/sources, L is the
  mixing matrix, M is unity), "factor analysis" (Z contains factors, L the loadings,
  M is unity, U the noise, Psi the noise covariance, lapla is a
  variational parameter for non-Gaussian factors, avini and ini are the
  information the factors convey about the observations).      
  
 
  }

  \details{

  \describe{
    \item{\bold{Plot}}{Plots all haplotype clusters of a haplotype
    cluster list, where the SNVs within
    the cluster are shown. Minor alleles of tagSNVs are marked by a
    particular color. Minor
    alleles of other SNVs (private, common, from other clusters, etc) are
    marked by a different color.  The model from \pkg{fabia} is also
    shown. 
    Asks for the next plot by pushing a key.
}



\item{\bold{Summary}}{Prints the number of haplotype clusters in the
  list and some statistics if available.}

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

data(hapRes)

mergedHaploClusterList <-
   hapRes$mergedHaploClusterList

# Summary method
summary(mergedHaploClusterList)

# Accessors
lengthList(mergedHaploClusterList)

statistics(mergedHaploClusterList)

summary(haploClusters(mergedHaploClusterList))


# Subsets

summary(mergedHaploClusterList[[1]])

summary(mergedHaploClusterList[1])

mergedHaploClusterList[[2]] <-
   mergedHaploClusterList[[1]]
mergedHaploClusterList[[3]] <-
   hapRes$mergedHaploClusterList2[[1]]

summary(mergedHaploClusterList)

mergedHaploClusterList[c(3,4)] <-
   mergedHaploClusterList[c(1,2)]

summary(mergedHaploClusterList)



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
