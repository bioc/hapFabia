\name{hapRes}
\docType{data}
\alias{hapRes}
\title{Example result of \code{hapFabia} }
\description{

  Results of a \code{hapFabia} call.
  List containing
 \enumerate{
 \item \code{mergedIBDsegmentList}: an object of the class
   \code{IBDsegmentList} that contains the extracted IBD segments
   that were extracted from two histograms with different offset.
 \item \code{res}: the result of FABIA.
 \item \code{sPF}: samples per loading of this FABIA result.
 \item \code{annot}: annotation for the genotype data.
 \item \code{IBDsegmentList1}: an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the first histogram.
 \item \code{IBDsegmentList2}: an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the second histogram.
 \item \code{mergedIBDsegmentList1}: an object of the class
   \code{IBDsegmentList} that contains the merged result of the first IBD segment
   extraction (redundancies removed).
 \item \code{mergedIBDsegmentList2}: an object of the class
   \code{IBDsegmentList} that contains the merged result of the second IBD segment
   extraction (redundancies removed).

}
  
}
\usage{hapRes}
\format{  List containing
 \enumerate{
 \item \code{mergedIBDsegmentList}: an object of the class
   \code{IBDsegmentList} that contains the extracted IBD segments
   that were extracted from two histograms with different offset.
 \item \code{res}: the result of FABIA.
 \item \code{sPF}: samples per loading of this FABIA result.
 \item \code{annot}: annotation for the genotype data.
 \item \code{IBDsegmentList1}: an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the first histogram.
 \item \code{IBDsegmentList2}: an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the second histogram.
 \item \code{mergedIBDsegmentList1}: an object of the class
   \code{IBDsegmentList} that contains the merged result of the first IBD segment
   extraction (redundancies removed).
 \item \code{mergedIBDsegmentList2}: an object of the class
   \code{IBDsegmentList} that contains the merged result of the second IBD segment
   extraction (redundancies removed).
}
}
\source{result of a call of \code{hapFabia} }
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
\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{datasets}


