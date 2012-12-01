\name{IBDsegmentList2excel}

\alias{IBDsegmentList2excel}
\alias{IBDsegmentList2excel,IBDsegmentList-method}
\alias{IBDsegmentList2excel,IBDsegmentList,character-method}

\title{Store an IBD segment list in EXCEL / csv format}

\description{

  \code{IBDsegmentList2excel}: \R implementation of \code{IBDsegmentList2excel}.


  IBD segment list is stored in a file in EXCEL format, more precise
  in comma separated format (.csv).

}
\usage{

\S4method{IBDsegmentList2excel}{IBDsegmentList,character}(IBDsegmentList,filename)


}
\arguments{
  \item{IBDsegmentList}{list of IBD segments given as an object
    of the class \code{IBDsegmentList}.}
  \item{filename}{name of the file where the IBD segment list is
    stored in EXCEL format.}

}
\details{


  IBD segment list is stored
  in comma separate format (\code{.csv}) which can readily be read by EXCEL.

  The EXCEL (\code{.csv}) file contains following columns:

\enumerate{


  \item{} \code{ID}: number of the IBD segment in the current extraction.
  \item{} \code{bicluster_id}: ID of the bicluster the IBD segment was found in.
  \item{} \code{chromosome}: the chromosome.
  \item{} \code{IBDsegmentPos}: genomic location of the IBD segment.
  \item{} \code{IBDsegmentLength}: length of the IBD segment.
  \item{} \code{numberIndividuals}: number of samples belonging to the IBD segment.
  \item{} \code{numbertagSNVs}: number tagSNVs marking the IBD segment.
  \item{} \code{individuals}: IDs of individuals or chromosomes belonging to the IBD segment.
  \item{} \code{tagSNVs}: IDs of SNVs that mark the IBD segment (tagSNVs).
  \item{} \code{populationIndividuals}: the population each individual belongs to.
  \item{} \code{idIndividuals}: IDs of the individuals or chromosomes.
  \item{} \code{labelIndividuals}: label of the individuals.
  \item{} \code{platformIndividuals}: for each individual the technology/platform that was used to genotype it.
  \item{} \code{coreClusterIndividuals}: IDs of individuals that constitute the core
  of the IBD segment.
  \item{} \code{tagSNVPositions}: physical positions of the tagSNVs on the
  chromosome in base pairs.
  \item{} \code{tagSNVAlleles}: alleles of the tagSNVs in the form \code{Ref:Alt} where Ref denotes reference allele and Alt the alternative allele.
  \item{} \code{tagSNVNames}: name of the tagSNVs according to a given annotation.
  \item{} \code{tagSNVFreq}: frequency of tagSNVs in the whole data set.
  \item{} \code{tagSNVGroupFreq}: frequency of tagSNVs in the population that is considered.
  \item{} \code{tagSNVChange}: if the minor allele was more frequent than the
  major, then both were switched. Switching is marked by a 1 and
  otherwise it is 0.
  \item{} \code{tagSNVsPerIndividual}: for each sample: tagSNVs are counted for which the sample has the minor allele.
  \item{} \code{individualPerTagSNV}: for each tagSNV: samples are counted for which the SNV has its minor allele.
\item{} \code{tagSNVAnno}: the functional annotation of tagSNVs for each tagSNV: like stop-loss,
  stop-gain, non-synonymous, synonymous, promoter, exonic, intronic, intergenic, etc.

}

  Implementation in \R.


}

\value{
writes to comma separated \code{.csv} file  
}

\examples{
old_dir <- getwd()
setwd(tempdir())

data(hapRes)
mergedIBDsegmentList <- hapRes$mergedIBDsegmentList
IBDsegmentList2excel(IBDsegmentList=mergedIBDsegmentList,
   filename="testResult.csv")

setwd(old_dir)

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

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}

\keyword{file}

\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare IBD segments,short IBD segments}
