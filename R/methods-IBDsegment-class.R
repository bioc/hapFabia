### ---------------------------
### IBDsegment class methods
### ---------------------------


##
## Constructor
##

IBDsegment <- function(ID=0,bicluster_id=0,chromosome="",IBDsegmentPos=0,IBDsegmentLength=0,numberIndividuals=0,numbertagSNVs=0,individuals=as.vector(0),tagSNVs=as.vector(0),populationIndividuals=as.vector(""),idIndividuals=as.vector(0),labelIndividuals=as.vector(""),platformIndividuals=as.vector(""),coreClusterIndividuals=as.vector(0),tagSNVPositions=as.vector(0),tagSNVAlleles=as.vector(""),tagSNVNames=as.vector(""),tagSNVFreq=as.vector(0),tagSNVGroupFreq=as.vector(0),tagSNVChange=as.vector(0),tagSNVsPerIndividual=as.vector(0),individualPerTagSNV=as.vector(0),tagSNVAnno=as.vector(c(list("")))) {

   new("IBDsegment",ID=ID,bicluster_id=bicluster_id,chromosome=chromosome,IBDsegmentPos=IBDsegmentPos,IBDsegmentLength=IBDsegmentLength,numberIndividuals=numberIndividuals,numbertagSNVs=numbertagSNVs,individuals=individuals,tagSNVs=tagSNVs,populationIndividuals=populationIndividuals,idIndividuals=idIndividuals,labelIndividuals=labelIndividuals,platformIndividuals=platformIndividuals,coreClusterIndividuals=coreClusterIndividuals,tagSNVPositions=tagSNVPositions,tagSNVAlleles=tagSNVAlleles,tagSNVNames=tagSNVNames,tagSNVFreq=tagSNVFreq,tagSNVGroupFreq=tagSNVGroupFreq,tagSNVChange=tagSNVChange,tagSNVsPerIndividual=tagSNVsPerIndividual,individualPerTagSNV=individualPerTagSNV,tagSNVAnno=tagSNVAnno)
}


##
## Getters and setters
##


setMethod("ID", "IBDsegment",
    function(x)
    {
      slot(x, "ID")
    }
)


setReplaceMethod("ID", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "ID") <- value
       x
    }
)


setMethod("bicluster_id", "IBDsegment",
    function(x)
    {
      slot(x, "bicluster_id")
    }
)


setReplaceMethod("bicluster_id", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "bicluster_id") <- value
       x
    }
)


setMethod("chromosome", "IBDsegment",
    function(x)
    {
      slot(x, "chromosome")
    }
)


setReplaceMethod("chromosome", c("IBDsegment", "character"),
    function(x, value)
    {
       slot(x, "chromosome") <- value
       x
    }
)


setMethod("IBDsegmentPos", "IBDsegment",
    function(x)
    {
      slot(x, "IBDsegmentPos")
    }
)


setReplaceMethod("IBDsegmentPos", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "IBDsegmentPos") <- value
       x
    }
)


setMethod("IBDsegmentLength", "IBDsegment",
    function(x)
    {
      slot(x, "IBDsegmentLength")
    }
)


setReplaceMethod("IBDsegmentLength", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "IBDsegmentLength") <- value
       x
    }
)


setMethod("numberIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "numberIndividuals")
    }
)


setReplaceMethod("numberIndividuals", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "numberIndividuals") <- value
       x
    }
)


setMethod("numbertagSNVs", "IBDsegment",
    function(x)
    {
      slot(x, "numbertagSNVs")
    }
)


setReplaceMethod("numbertagSNVs", c("IBDsegment", "numeric"),
    function(x, value)
    {
       slot(x, "numbertagSNVs") <- value
       x
    }
)


setMethod("individuals", "IBDsegment",
    function(x)
    {
      slot(x, "individuals")
    }
)


setReplaceMethod("individuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "individuals") <- value
       x
    }
)


setMethod("tagSNVs", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVs")
    }
)


setReplaceMethod("tagSNVs", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVs") <- value
       x
    }
)


setMethod("populationIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "populationIndividuals")
    }
)


setReplaceMethod("populationIndividuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "populationIndividuals") <- value
       x
    }
)


setMethod("idIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "idIndividuals")
    }
)


setReplaceMethod("idIndividuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "idIndividuals") <- value
       x
    }
)


setMethod("labelIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "labelIndividuals")
    }
)


setReplaceMethod("labelIndividuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "labelIndividuals") <- value
       x
    }
)


setMethod("platformIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "platformIndividuals")
    }
)


setReplaceMethod("platformIndividuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "platformIndividuals") <- value
       x
    }
)


setMethod("coreClusterIndividuals", "IBDsegment",
    function(x)
    {
      slot(x, "coreClusterIndividuals")
    }
)

setReplaceMethod("coreClusterIndividuals", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "coreClusterIndividuals") <- value
       x
    }
)


setMethod("tagSNVPositions", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVPositions")
    }
)


setReplaceMethod("tagSNVPositions", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVPositions") <- value
       x
    }
)


setMethod("tagSNVAlleles", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVAlleles")
    }
)


setReplaceMethod("tagSNVAlleles", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVAlleles") <- value
       x
    }
)


setMethod("tagSNVNames", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVNames")
    }
)


setReplaceMethod("tagSNVNames", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVNames") <- value
       x
    }
)


setMethod("tagSNVFreq", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVFreq")
    }
)


setReplaceMethod("tagSNVFreq", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVFreq") <- value
       x
    }
)


setMethod("tagSNVGroupFreq", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVGroupFreq")
    }
)


setReplaceMethod("tagSNVGroupFreq", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVGroupFreq") <- value
       x
    }
)



setMethod("tagSNVChange", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVChange")
    }
)


setReplaceMethod("tagSNVChange", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVChange") <- value
       x
    }
)


setMethod("tagSNVsPerIndividual", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVsPerIndividual")
    }
)


setReplaceMethod("tagSNVsPerIndividual", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVsPerIndividual") <- value
       x
    }
)


setMethod("individualPerTagSNV", "IBDsegment",
    function(x)
    {
      slot(x, "individualPerTagSNV")
    }
)


setReplaceMethod("individualPerTagSNV", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "individualPerTagSNV") <- value
       x
    }
)


setMethod("tagSNVAnno", "IBDsegment",
    function(x)
    {
      slot(x, "tagSNVAnno")
    }
)


setReplaceMethod("tagSNVAnno", c("IBDsegment", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVAnno") <- value
       x
    }
)





##
## Summary
##



setMethod("summary", "IBDsegment",
function(object, ...)
{
    cat("\nAn object of class",class(object))

    cat("\nIBD segment ID: ",ID(object))
    cat("\nFrom bicluster: ",bicluster_id(object))
    cat("\nChromosome: ",chromosome(object))
    cat("\nPosition: ",prettyNum(IBDsegmentPos(object), big.mark = ","))
    cat("\nLength SNVs: ",round(IBDsegmentLength(object)))
    cat("\nLength: ",round((max(tagSNVPositions(object))-min(tagSNVPositions(object)))/1000) , "kbp")
    cat("\nNumber of individuals/chromosomes: ",numberIndividuals(object))
    cat("\nNumber of tagSNVs: ",numbertagSNVs(object), "\n")

})




##
## Plot
##

setMethod("plot",signature(x="IBDsegment", y="missing"),
function(x,filename, ...) {


    if (missing(x)) {
        stop("IBD segment 'x' is missing. Stopped.")
    }

    if (missing(filename)) {
        stop("File name 'filename' for loading genotype data is missing. Stopped.")
    }




individ <- individuals(x)
tagSNV <- tagSNVs(x)
tagSNV <- as.integer(sort.int(as.integer(unique(tagSNV))))
tagSNVPositions <- tagSNVPositions(x)
chrom <- chromosome(x)

labels_ALL <- labelIndividuals(x)

Lout <- readSamplesSpfabia(X=filename,samples=individ,lowerB=0,upperB=1000.0)


tagSNVL <- list(tagSNV)
labelsK <- c("model L")

# here other annotations can be included
# tagSNVL <- list(tagSNV,tagSNVAnnot1,tagSNVAnnot2)
# labelsK <- c("model L","Annot1",Annot2")

plotIBDsegment(Lout=Lout,tagSNV=tagSNVL,physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),chrom=chrom,count=0,labelsNA=labels_ALL,labelsNA1=labelsK,...)


}
)


setMethod("plotLarger",signature(x="IBDsegment", filename="character",fact="numeric",addSamp="ANY"),
function(x,filename,fact=1.0,addSamp=c(), ...) {


    if (missing(x)) {
        stop("IBD segment 'x' is missing. Stopped.")
    }

    if (missing(filename)) {
        stop("File name 'filename' for loading genotype data is missing. Stopped.")
    }




individ <- individuals(x)
intss <- intersect(individ,addSamp)
addSamp <- setdiff(addSamp,intss)
samp <- c(addSamp,individ)
l1 <- length(individ)
l2 <- length(addSamp)
indik <- c(rep("A",l2),rep("B",l1))

labels1 <- labelIndividuals(x)
labels2 <- as.character(addSamp)
labels12 <- c(labels2,labels1)
labels_ALL <- paste(labels12,indik,sep="_")

sortS <- sort.int(as.integer(unique(samp)),index.return = TRUE)
individ <- samp[sortS$ix]
labels_ALL <- labels_ALL[sortS$ix]



tagSNV <- tagSNVs(x)
tagSNV <- as.integer(sort.int(as.integer(unique(tagSNV))))
tagSNVPositions <- tagSNVPositions(x)



ltagSNV <- length(tagSNV)
llA <- tagSNV[ltagSNV]-tagSNV[1]
llAd <- round((llA*fact-llA)/2)
tagSNVEnd <- tagSNV[ltagSNV]+llAd
tagSNVBeg <- tagSNV[1]-llAd
tagSNV <- c(tagSNVBeg,tagSNV,tagSNVEnd)
tagSNVPositions <- c(tagSNVPositions[1],tagSNVPositions,tagSNVPositions[ltagSNV])


chrom <- chromosome(x)


Lout <- readSamplesSpfabia(X=filename,samples=individ,lowerB=0,upperB=1000.0)


tagSNVL <- list(tagSNV)
labelsK <- c("model L")

# here other annotations can be included
# tagSNVL <- list(tagSNV,tagSNVAnnot1,tagSNVAnnot2)
# labelsK <- c("model L","Annot1",Annot2")

plotIBDsegment(Lout=Lout,tagSNV=tagSNVL,physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),chrom=chrom,count=0,labelsNA=labels_ALL,labelsNA1=labelsK,...)


}
)


