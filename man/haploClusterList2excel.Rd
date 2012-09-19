\name{haploClusterList2excel}

\alias{haploClusterList2excel}
\alias{haploClusterList2excel,HaploClusterList-method}
\alias{haploClusterList2excel,HaploClusterList,character-method}

\title{Store a haplotype cluster list in EXCEL / csv format}

\description{

  \code{haploClusterList2excel}: \R implementation of \code{haploClusterList2excel}.


  Haplotype cluster list is stored in a file in EXCEL format, more precise
  in comma separated format (.csv).

}
\usage{

\S4method{haploClusterList2excel}{HaploClusterList,character}(haploClusterList,filename)


}
\arguments{
  \item{haploClusterList}{list of haplotype clusters given as an object
    of the class \code{HaploClusterList}.}
  \item{filename}{name of the file where the haplotype cluster list is
    stored in EXCEL format.}

}
\details{


  Haplotype cluster list is stored
  in comma separate format (\code{.csv}) which can readily be read by EXCEL.

  The EXCEL (\code{.csv}) file contains following columns:

\enumerate{


  \item{} \code{ID}: number of the haplotype cluster in the current extraction.
  \item{} \code{bicluster_id}: ID of the bicluster the haplotype cluster was found in.
  \item{} \code{chromosome}: the chromosome.
  \item{} \code{haploClusterPos}: genomic location of the haplotype cluster.
  \item{} \code{haploClusterLength}: length of the haplotype cluster.
  \item{} \code{numberIndividuals}: number of samples belonging to the haplotype cluster.
  \item{} \code{numbertagSNVs}: number tagSNVs marking the haplotype cluster.
  \item{} \code{individuals}: IDs of individuals or chromosomes belonging to the haplotype cluster.
  \item{} \code{tagSNVs}: IDs of SNVs that mark the haplotype cluster (tagSNVs).
  \item{} \code{populationIndividuals}: the population each individual belongs to.
  \item{} \code{idIndividuals}: IDs of the individuals or chromosomes.
  \item{} \code{labelIndividuals}: label of the individuals.
  \item{} \code{platformIndividuals}: for each individual the technology/platform that was used to genotype it.
  \item{} \code{coreClusterIndividuals}: IDs of individuals that constitute the core
  of the haplotype cluster.
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
mergedHaploClusterList <- hapRes$mergedHaploClusterList
haploClusterList2excel(haploClusterList=mergedHaploClusterList,
   filename="testResult.csv")

setwd(old_dir)

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

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}

\keyword{file}

\concept{genetics,haplotype,identity by descent,bicluster,next generation sequencing,genotype,single nucleotide polymorphism,single nucleotide variation,rare variants,rare SNPs, rare SNVs,rare haplotype clusters,short haplotype clusters}
