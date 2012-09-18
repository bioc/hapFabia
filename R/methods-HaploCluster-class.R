### ---------------------------
### HaploCluster class methods
### ---------------------------


##
## Constructor
##

HaploCluster <- function(ID=0,bicluster_id=0,chromosome="",haploClusterPos=0,haploClusterLength=0,numberIndividuals=0,numbertagSNVs=0,individuals=as.vector(0),tagSNVs=as.vector(0),populationIndividuals=as.vector(""),idIndividuals=as.vector(0),labelIndividuals=as.vector(""),platformIndividuals=as.vector(""),coreClusterIndividuals=as.vector(0),tagSNVPositions=as.vector(0),tagSNVAlleles=as.vector(""),tagSNVNames=as.vector(""),tagSNVFreq=as.vector(0),tagSNVGroupFreq=as.vector(0),tagSNVChange=as.vector(0),tagSNVsPerIndividual=as.vector(0),individualPerTagSNV=as.vector(0),tagSNVAnno=as.vector(c(list("")))) {

   new("HaploCluster",ID=ID,bicluster_id=bicluster_id,chromosome=chromosome,haploClusterPos=haploClusterPos,haploClusterLength=haploClusterLength,numberIndividuals=numberIndividuals,numbertagSNVs=numbertagSNVs,individuals=individuals,tagSNVs=tagSNVs,populationIndividuals=populationIndividuals,idIndividuals=idIndividuals,labelIndividuals=labelIndividuals,platformIndividuals=platformIndividuals,coreClusterIndividuals=coreClusterIndividuals,tagSNVPositions=tagSNVPositions,tagSNVAlleles=tagSNVAlleles,tagSNVNames=tagSNVNames,tagSNVFreq=tagSNVFreq,tagSNVGroupFreq=tagSNVGroupFreq,tagSNVChange=tagSNVChange,tagSNVsPerIndividual=tagSNVsPerIndividual,individualPerTagSNV=individualPerTagSNV,tagSNVAnno=tagSNVAnno)
}


##
## Getters and setters
##


setMethod("ID", "HaploCluster",
    function(x)
    {
      slot(x, "ID")
    }
)


setReplaceMethod("ID", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "ID") <- value
       x
    }
)


setMethod("bicluster_id", "HaploCluster",
    function(x)
    {
      slot(x, "bicluster_id")
    }
)


setReplaceMethod("bicluster_id", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "bicluster_id") <- value
       x
    }
)


setMethod("chromosome", "HaploCluster",
    function(x)
    {
      slot(x, "chromosome")
    }
)


setReplaceMethod("chromosome", c("HaploCluster", "character"),
    function(x, value)
    {
       slot(x, "chromosome") <- value
       x
    }
)


setMethod("haploClusterPos", "HaploCluster",
    function(x)
    {
      slot(x, "haploClusterPos")
    }
)


setReplaceMethod("haploClusterPos", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "haploClusterPos") <- value
       x
    }
)


setMethod("haploClusterLength", "HaploCluster",
    function(x)
    {
      slot(x, "haploClusterLength")
    }
)


setReplaceMethod("haploClusterLength", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "haploClusterLength") <- value
       x
    }
)


setMethod("numberIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "numberIndividuals")
    }
)


setReplaceMethod("numberIndividuals", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "numberIndividuals") <- value
       x
    }
)


setMethod("numbertagSNVs", "HaploCluster",
    function(x)
    {
      slot(x, "numbertagSNVs")
    }
)


setReplaceMethod("numbertagSNVs", c("HaploCluster", "numeric"),
    function(x, value)
    {
       slot(x, "numbertagSNVs") <- value
       x
    }
)


setMethod("individuals", "HaploCluster",
    function(x)
    {
      slot(x, "individuals")
    }
)


setReplaceMethod("individuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "individuals") <- value
       x
    }
)


setMethod("tagSNVs", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVs")
    }
)


setReplaceMethod("tagSNVs", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVs") <- value
       x
    }
)


setMethod("populationIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "populationIndividuals")
    }
)


setReplaceMethod("populationIndividuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "populationIndividuals") <- value
       x
    }
)


setMethod("idIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "idIndividuals")
    }
)


setReplaceMethod("idIndividuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "idIndividuals") <- value
       x
    }
)


setMethod("labelIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "labelIndividuals")
    }
)


setReplaceMethod("labelIndividuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "labelIndividuals") <- value
       x
    }
)


setMethod("platformIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "platformIndividuals")
    }
)


setReplaceMethod("platformIndividuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "platformIndividuals") <- value
       x
    }
)


setMethod("coreClusterIndividuals", "HaploCluster",
    function(x)
    {
      slot(x, "coreClusterIndividuals")
    }
)

setReplaceMethod("coreClusterIndividuals", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "coreClusterIndividuals") <- value
       x
    }
)


setMethod("tagSNVPositions", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVPositions")
    }
)


setReplaceMethod("tagSNVPositions", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVPositions") <- value
       x
    }
)


setMethod("tagSNVAlleles", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVAlleles")
    }
)


setReplaceMethod("tagSNVAlleles", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVAlleles") <- value
       x
    }
)


setMethod("tagSNVNames", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVNames")
    }
)


setReplaceMethod("tagSNVNames", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVNames") <- value
       x
    }
)


setMethod("tagSNVFreq", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVFreq")
    }
)


setReplaceMethod("tagSNVFreq", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVFreq") <- value
       x
    }
)


setMethod("tagSNVGroupFreq", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVGroupFreq")
    }
)


setReplaceMethod("tagSNVGroupFreq", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVGroupFreq") <- value
       x
    }
)



setMethod("tagSNVChange", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVChange")
    }
)


setReplaceMethod("tagSNVChange", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVChange") <- value
       x
    }
)


setMethod("tagSNVsPerIndividual", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVsPerIndividual")
    }
)


setReplaceMethod("tagSNVsPerIndividual", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVsPerIndividual") <- value
       x
    }
)


setMethod("individualPerTagSNV", "HaploCluster",
    function(x)
    {
      slot(x, "individualPerTagSNV")
    }
)


setReplaceMethod("individualPerTagSNV", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "individualPerTagSNV") <- value
       x
    }
)


setMethod("tagSNVAnno", "HaploCluster",
    function(x)
    {
      slot(x, "tagSNVAnno")
    }
)


setReplaceMethod("tagSNVAnno", c("HaploCluster", "vector"),
    function(x, value)
    {
       slot(x, "tagSNVAnno") <- value
       x
    }
)





##
## Summary
##



setMethod("summary", "HaploCluster",
function(object, ...)
{
    cat("\nAn object of class",class(object))

    cat("\nHaplotype cluster ID: ",ID(object))
    cat("\nFrom bicluster: ",bicluster_id(object))
    cat("\nChromosome: ",chromosome(object))
    cat("\nPosition: ",prettyNum(haploClusterPos(object), big.mark = ","))
    cat("\nLength SNVs: ",round(haploClusterLength(object)))
    cat("\nLength: ",round((max(tagSNVPositions(object))-min(tagSNVPositions(object)))/1000) , "kbp")
    cat("\nNumber of individuals/chromosomes: ",numberIndividuals(object))
    cat("\nNumber of tagSNVs: ",numbertagSNVs(object), "\n")

})




##
## Plot
##

setMethod("plot",signature(x="HaploCluster", y="missing"),
function(x,filename, ...) {

    require(fabia)

    if (missing(x)) {
        stop("Haplotype cluster 'x' is missing. Stopped.")
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

plotHaplotypeCluster(Lout=Lout,tagSNV=tagSNVL,physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),chrom=chrom,count=0,labelsNA=labels_ALL,labelsNA1=labelsK)


}
)

