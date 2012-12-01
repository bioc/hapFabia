\name{hapFabia}

\alias{hapFabia}

\title{IBD segment extraction by FABIA}

\description{

  \code{hapFabia}: \R implementation of the \emph{hapFabia} method.

  \code{hapFabia} extracts \bold{short IBD segments} tagged by rare variants
  from phased or unphased genotypes. hapFabia is designed for rare variants
  in very large sequencing data.
  The method is based on FABIA biclustering and utilizes the package
  \pkg{fabia}.

}
\usage{

hapFabia(fileName,prefixPath="",sparseMatrixPostfix="_mat",
   annotPostfix="_annot.txt",individualsPostfix="_individuals.txt",
   labelsA=NULL,pRange="",individuals=0,lowerBP=0,upperBP=0.05,
   p=10,iter=40,quant=0.01,eps=1e-5,alpha=0.03,cyc=50,non_negative=1,
   write_file=0,norm=0,lap=100.0,IBDsegmentLength=50,Lt = 0.1,
   Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.03,
   haplotypes=FALSE,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,
   simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)
}


\arguments{

 \item{fileName}{name of the file that contains the genotype matrix in sparse format.}
 \item{prefixPath}{path to the genotype file.}
 \item{sparseMatrixPostfix}{postfix string for the sparse matrix.}
 \item{annotPostfix}{postfix string for the SNV annotation file.}
 \item{individualsPostfix}{postfix string for the file containing the
 names of the individuals.}
 \item{labelsA}{annotation of the individuals as matrix individuals x 4 
(individual ID,
 subpopulation, population, genotyping platform).}
 \item{pRange}{indicates which DNA interval is processed.}
 \item{individuals}{vector of individuals that are included into the analysis; default = 0 (all individuals).}
 \item{lowerBP}{lower bound on minor allele frequencies (MAF); however at least two occurrences are required to remove private SNVs.}
 \item{upperBP}{upper bound  on minor allele frequencies (MAF) to extract rare variants.}
 \item{p}{number of biclusters per fabia iteration.}
 \item{iter}{number of fabia iterations.}
 \item{quant}{percentage of fabia loadings L that are removed after each iteration.}
 \item{eps}{lower bound on fabia variational parameter lapla; default 1e-5.}
 \item{alpha}{fabia sparseness of the loadings; default = 0.03.}
 \item{cyc}{number of cycles per fabia iteration; default 50.}
 \item{non_negative}{non-negative fabia factors and loadings if non_negative = 1; default = 1 (yes).}
 \item{write_file}{fabia results are written to files (L in sparse format), default = 0 (not written).}
 \item{norm}{fabia data normalization; default 1 (no normalization).}
 \item{lap}{minimal value of fabia's variational parameter; default 100.0.}
 \item{IBDsegmentLength}{expected typical IBD segment length in kbp.}
 \item{Lt}{percentage of largest fabia L values to consider for IBD segment extraction.}
 \item{Zt}{percentage of largest fabia Z values to consider for IBD segment extraction.}
 \item{thresCount}{p-value of random histogram hit; default 1e-5.}
 \item{mintagSNVsFactor}{percentage of the histogram tagSNVs count
   threshold (\code{mintagSNVs} in \code{extractIBDsegments});
   used to define minimal overlap of individual intervals in an IBD segment; default 3/4.}
 \item{pMAF}{averaged and corrected (for non-uniform distributions) minor allele frequency.}
\item{haplotypes}{\code{haplotypes = TRUE} indicates phased genotypes
 that is two chromosomes per individual otherwise unphased genotypes.}
\item{cut}{cutoff for merging IBD segments after a hierarchical
 clustering; default 0.8.}
\item{procMinIndivids}{Percentage of cluster individuals a tagSNV must
 tag to be considered as tagSNV for the IBD segment.}
\item{thresPrune}{Threshold for pruning border tagSNVs based on an
 exponential distribution where border tagSNVs with large distances to
 the next tagSNV are pruned.}
\item{simv}{Similarity measure for merging clusters: \code{"minD"} (percentage of smaller explained by larger set),
   \code{"jaccard"} (Jaccard index), \code{"dice"} (Dice index), or
   \code{"maxD"}; default \code{"minD"}.}
\item{minTagSNVs}{Minimum matching tagSNVs for cluster similarity;
 otherwise the similarity is set to zero.}
\item{minIndivid}{Minimum matching individuals for cluster similarity;
 otherwise the similarity is set to zero.}
 \item{avSNVsDist}{average distance between SNVs in base pairs - used
 together with \code{IBDsegmentLength} to compute the number of SNVs
 in the histogram bins; default=100.}
\item{SNVclusterLength}{if \code{IBDsegmentLength=0} then the number
 of SNVs in the histogram bins can be given directly; default 100.}
}
\details{

  This function uses a genotype matrix in sparse matrix format
  and extracts IBD segments by
 biclustering.
 First, it performs Fabia biclustering and then extracts local
 accumulations of loadings.
 Then it disentangles IBD segments and prunes off spurious correlated SNVs.
 Finally, it merges similar IBD segments to account for larger IBD
  segments that were broken during the analysis.

 Annotation file \code{..._annot.txt} for SNVs:
\enumerate{
 \item  first line: number individuals;
 \item  second line: number SNVs;
 \item  for each SNV a line containing following field that are blank separated:
 "chromosome", "physical position", "snvNames", "snvMajor", "snvMinor",
 "quality", "pass", "info of vcf file", "fields in vcf file",
 "frequency", "0/1: 1 is changed if major allele is actually minor
 allele".
 }
 
 \code{labelsA} is a matrix ("number individuals" x 4),
 where for each individual following characteristics are given:
 \enumerate{
   \item id;
   \item subPopulation;
   \item population;
   \item platform.
 }
 
 The probability of observing \eqn{k}{k} or more
 correlated SNVs in a histogram bin is computed.
 The minimal \eqn{k}{k} which pushes the
 probability below the threshold \code{thresCount} is used to find
 accumulation of correlated SNVs that are assumed to belong to a
 IBD segment.

 Let \eqn{p}{p} be the probability of a random
 minor allele match between \eqn{t}{t} individuals.
The probability of observing \eqn{k}{k} or more matches for
 \eqn{n}{n} SNVs in
 a histogram bin is given by one minus the binomial distribution
 \eqn{F(k;n,p)}{F(k;n,p)}:

 \deqn{1 \ - \ F(k-1;n,p) \ = \ \Pr(K \geq k) \ = \
   \sum_{i=k}^\infty \ {n\choose i} \ p^i \ (1-p)^{n-i}}{1 - F(k-1;n,p)
   =  Pr(K >= k) = SUM(i = k to infinity; (n over i)  power(p,i) power(1-p,n-i))}
If \eqn{q}{q} is the minor allele frequency (MAF) for one SNV,
the probability \eqn{p}{p} of observing the minor allele of this SNV in
all \eqn{t}{t} individuals is \eqn{p=q^t}{p=power(q,t)}.
The value \eqn{q}{q} is given by the parameter \code{pMAF}.

 Implementation in \R.

}
\value{
List containing
 \item{mergedIBDsegmentList}{an object of the class
   \code{IBDsegmentList} that contains the extracted IBD segments
   that were extracted from two histograms with different offset.}
 \item{res}{the result of FABIA.}
 \item{sPF}{individuals per loading of this FABIA result.}
 \item{annot}{annotation for the genotype data.}
 \item{IBDsegmentList1}{an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the first histogram.}
 \item{IBDsegmentList2}{an object of the class
   \code{IBDsegmentList} that contains the result of IBD segment extraction from the second histogram.}
 \item{mergedIBDsegmentList1}{an object of the class
   \code{IBDsegmentList} that contains the merged result of the first IBD segment
   extraction (redundancies removed).}
 \item{mergedIBDsegmentList2}{an object of the class
   \code{IBDsegmentList} that contains the merged result of the second IBD segment
   extraction (redundancies removed).}

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
   sep = " ",quote = FALSE,row.names = FALSE,
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


hapRes <- hapFabia(fileName="dataSim1fabia",prefixPath="",
   sparseMatrixPostfix="_mat",
   annotPostfix="_annot.txt",individualsPostfix="_individuals.txt",
   labelsA=NULL,pRange="",individuals=0,lowerBP=0,upperBP=0.15,
   p=10,iter=1,quant=0.01,eps=1e-5,alpha=0.03,cyc=50,non_negative=1,
   write_file=0,norm=0,lap=100.0,IBDsegmentLength=10,Lt = 0.1,
   Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.1,
   haplotypes=FALSE,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,
   simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)

summary(hapRes$mergedIBDsegmentList)

plot(hapRes$mergedIBDsegmentList[[1]],filename="dataSim1fabia_mat")


### Another Example
simulateIBDsegmentsFabia(fileprefix="dataSim",
   minruns=2,maxruns=2,snvs=1000,individualsN=100,
   avDistSnvs=100,avDistMinor=10,noImplanted=1,
   implanted=10,length=50,minors=30,mismatches=0,
   mismatchImplanted=0.5,overlap=50)

hapRes <- hapFabia(fileName="dataSim2fabia",prefixPath="",
   sparseMatrixPostfix="_mat",
   annotPostfix="_annot.txt",individualsPostfix="_individuals.txt",
   labelsA=NULL,pRange="",individuals=0,lowerBP=0,upperBP=0.15,
   p=10,iter=1,quant=0.01,eps=1e-5,alpha=0.03,cyc=50,non_negative=1,
   write_file=0,norm=0,lap=100.0,IBDsegmentLength=10,Lt = 0.1,
   Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.1,
   haplotypes=FALSE,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,
   simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)

## Summary of the IBD segment list
summary(hapRes$mergedIBDsegmentList)

## Summary of the IBD segment 
summary(hapRes$mergedIBDsegmentList[[1]])


## Plot an IBD segment
plot(hapRes$mergedIBDsegmentList[[1]],filename="dataSim2fabia_mat")

\dontrun{
## It is interactive, thus dontrun!

## Plot an IBD segment list
plot(hapRes$mergedIBDsegmentList,filename="dataSim2fabia_mat")

}

setwd(old_dir)


}

\references{

S. Hochreiter et al.,
\sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
Bioinformatics 26(12):1520-1527, 2010.

}
\keyword{models,multivariate,cluster}
\concept{genetics,haplotype,identity by descent,bicluster,next
  generation sequencing,genotype,single nucleotide polymorphism,single
  nucleotide variation,rare variants,rare SNVs, rare
  IBD segments,short IBD segments, IBD}
