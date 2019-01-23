\name{IBDsegmentList-class}
\docType{class}

% Class
\alias{IBDsegmentList-class}
\alias{IBDsegmentList-method}

% Constructor-like function:
\alias{IBDsegmentList}
\alias{IBDsegmentList,ANY-method}
\alias{IBDsegmentList,list,numeric,list-method}

% Coercion methods:

% Accessor methods:

\alias{IBDsegments}
\alias{IBDsegments,IBDsegmentList-method}
\alias{IBDsegments<-}
\alias{IBDsegments<-,IBDsegmentList,list-method}


\alias{lengthList}
\alias{lengthList,IBDsegmentList-method}
\alias{lengthList<-}
\alias{lengthList<-,IBDsegmentList,numeric-method}


\alias{statistics}
\alias{statistics,IBDsegmentList-method}
\alias{statistics<-}
\alias{statistics<-,IBDsegmentList,list-method}


\alias{[[}
\alias{[[,IBDsegmentList-method}
\alias{[[,IBDsegmentList,numeric,missing-method}
\alias{[[<-}
\alias{[[<-,IBDsegmentList-method}
\alias{[[<-,IBDsegmentList,numeric,missing,IBDsegment-method}


\alias{[}
\alias{[,IBDsegmentList-method}
\alias{[,IBDsegmentList,numeric,missing-method}
\alias{[<-}
\alias{[<-,IBDsegmentList-method}
\alias{[<-,IBDsegmentList,numeric,missing,IBDsegmentList-method}



% summary method:
\alias{summary,IBDsegmentList-method}
% plot method:
\alias{plot,IBDsegmentList-method}
\alias{plot,IBDsegmentList,missing-method}



\title{IBDsegmentList instances and methods}

\description{\code{IBDsegmentList} is a class to store a list of
  IBD segments with its statistics.
  Lists can be merged or analyzed in subsequent steps.
}

\section{Slots}{
  Objects of class \code{IBDsegmentList} have the following slots:
  \describe{
        \item{\code{IBDsegments}:}{List of IBD segments.}
         \item{\code{lengthList}:}{Number of IBD segments in the list.}
         \item{\code{statistics}:}{Statistics of IBD segments like
	   average length, average number of individuals belonging to
	   an IBD segment, average number of tagSNVs of an IBD segment, etc.}

     }
  }

\section{Constructor}{
 Constructor of class IBDsegmentList.
 
 \code{IBDsegmentList(IBDsegments=list(),lengthList=0,statistics=list())}
   
}

\section{Accessors}{
  In the following \code{x} denotes an IBDsegmentList object.
   
 
  \describe{

    \item{}{
      \code{IBDsegments(x)}, \code{IBDsegments(x) <- value}:
      Returns or sets \code{IBDsegments}, where the return value and
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
      \code{value} are both an instance of the class \code{IBDsegment}. 
    }

    \item{}{
      \code{x[i]}, \code{x[i] <- value}:
      Returns or sets a sublist of the list \code{x},
      where the return value and
      \code{value} are both an instance of the class \code{IBDsegmentList}. 
    }


  }
}


\usage{


\S4method{plot}{IBDsegmentList}(x, ... )


\S4method{summary}{IBDsegmentList}(object, ...) 

}


\section{Signatures}{
  \describe{


\item{plot}{\code{signature(x = "IBDsegmentList", y = "missing")}}{
  Plotting all IBD segments of the list using an interactive command.}


\item{summary}{\code{signature(object = "IBDsegmentList")}}{ Summary
  of a list of IBD segments, where the number of clusters and a
  statistics are given. }


}
}

\arguments{
  \item{x}{object of the class \code{IBDsegment}.}
  \item{object}{object of the class \code{IBDsegment}.}
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
  IBDsegmentList objects are returned by \code{fabia}, \code{fabias}, \code{fabiap},
  \code{fabiasp},  \code{mfsc}, \code{nmfsc},
  \code{nmfdiv}, and \code{nmfeu}.} 

\section{Extension to store results of other methods}{
  
  The class \code{IBDsegmentList} may contain the result of different matrix factorization
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
    \item{\bold{Plot}}{Plots all IBD segments of
      an IBD segment list, where the SNVs within
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
    Asks for the next plot by pushing a key.
}





\item{\bold{Summary}}{Prints the number of IBD segments in the
  list and some statistics if available.}

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

data(hapRes)

mergedIBDsegmentList <-
   hapRes$mergedIBDsegmentList

# Summary method
summary(mergedIBDsegmentList)

# Accessors
lengthList(mergedIBDsegmentList)

statistics(mergedIBDsegmentList)

summary(IBDsegments(mergedIBDsegmentList))


# Subsets

summary(mergedIBDsegmentList[[1]])

summary(mergedIBDsegmentList[1])

mergedIBDsegmentList[[2]] <-
   mergedIBDsegmentList[[1]]
mergedIBDsegmentList[[3]] <-
   hapRes$mergedIBDsegmentList2[[1]]

summary(mergedIBDsegmentList)

# mergedIBDsegmentList[c(3,4)] <-
#    mergedIBDsegmentList[c(1,2)]
# summary(mergedIBDsegmentList)
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
