#
#
# Author: SEPP HOCHREITER
###############################################################################


## class Factorization

setGeneric("topLZ",
           signature = c("res","n","LZ","indices","p","w"),
           function(res,n=1,LZ="L",indices=TRUE,p=NULL,w=NULL)
           standardGeneric("topLZ"))


setGeneric("plotL",
           signature = c("res","n","p","w","type","intervv","off","t","cex"),
           function(res,n=1,p=NULL,w=NULL,type="points",intervv=500,off=0,t="p",cex=1)
           standardGeneric("plotL"))


setGeneric("histL",
           signature = c("res","n","p","w","intervv","off"),
           function(res,n=1,p=NULL,w=NULL,intervv=500,off=0)
           standardGeneric("histL"))


setGeneric("extractHaploClusters",
           signature = c("res","sPF","annot","chrom","labelsA","ps","psZ","inteA","thresA","mintagSNVs","off","procMinIndivids","thresPrune"),
           function(res,sPF,annot,chrom,labelsA,ps,psZ,inteA,thresA,mintagSNVs,off,procMinIndivids,thresPrune)
           standardGeneric("extractHaploClusters"))



## class HaploCluster



setGeneric("ID", signature = "x",
    function(x) standardGeneric("ID")
)

setGeneric("ID<-", signature = c("x", "value"),
    function(x, value) standardGeneric("ID<-")
)

setGeneric("bicluster_id", signature = "x",
    function(x) standardGeneric("bicluster_id")
)

setGeneric("bicluster_id<-", signature = c("x", "value"),
    function(x, value) standardGeneric("bicluster_id<-")
)

setGeneric("chromosome", signature = "x",
    function(x) standardGeneric("chromosome")
)

setGeneric("chromosome<-", signature = c("x", "value"),
    function(x, value) standardGeneric("chromosome<-")
)

setGeneric("haploClusterPos", signature = "x",
    function(x) standardGeneric("haploClusterPos")
)

setGeneric("haploClusterPos<-", signature = c("x", "value"),
    function(x, value) standardGeneric("haploClusterPos<-")
)

setGeneric("haploClusterLength", signature = "x",
    function(x) standardGeneric("haploClusterLength")
)

setGeneric("haploClusterLength<-", signature = c("x", "value"),
    function(x, value) standardGeneric("haploClusterLength<-")
)

setGeneric("numberIndividuals", signature = "x",
    function(x) standardGeneric("numberIndividuals")
)

setGeneric("numberIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("numberIndividuals<-")
)

setGeneric("numbertagSNVs", signature = "x",
    function(x) standardGeneric("numbertagSNVs")
)

setGeneric("numbertagSNVs<-", signature = c("x", "value"),
    function(x, value) standardGeneric("numbertagSNVs<-")
)

setGeneric("individuals", signature = "x",
    function(x) standardGeneric("individuals")
)

setGeneric("individuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("individuals<-")
)

setGeneric("tagSNVs", signature = "x",
    function(x) standardGeneric("tagSNVs")
)

setGeneric("tagSNVs<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVs<-")
)

setGeneric("populationIndividuals", signature = "x",
    function(x) standardGeneric("populationIndividuals")
)

setGeneric("populationIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("populationIndividuals<-")
)

setGeneric("idIndividuals", signature = "x",
    function(x) standardGeneric("idIndividuals")
)

setGeneric("idIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("idIndividuals<-")
)

setGeneric("labelIndividuals", signature = "x",
    function(x) standardGeneric("labelIndividuals")
)

setGeneric("labelIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("labelIndividuals<-")
)

setGeneric("platformIndividuals", signature = "x",
    function(x) standardGeneric("platformIndividuals")
)

setGeneric("platformIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("platformIndividuals<-")
)

setGeneric("coreClusterIndividuals", signature = "x",
    function(x) standardGeneric("coreClusterIndividuals")
)

setGeneric("coreClusterIndividuals<-", signature = c("x", "value"),
    function(x, value) standardGeneric("coreClusterIndividuals<-")
)

setGeneric("tagSNVPositions", signature = "x",
    function(x) standardGeneric("tagSNVPositions")
)

setGeneric("tagSNVPositions<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVPositions<-")
)

setGeneric("tagSNVAlleles", signature = "x",
    function(x) standardGeneric("tagSNVAlleles")
)

setGeneric("tagSNVAlleles<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVAlleles<-")
)


setGeneric("tagSNVNames", signature = "x",
    function(x) standardGeneric("tagSNVNames")
)

setGeneric("tagSNVNames<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVNames<-")
)


setGeneric("tagSNVFreq", signature = "x",
    function(x) standardGeneric("tagSNVFreq")
)

setGeneric("tagSNVFreq<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVFreq<-")
)


setGeneric("tagSNVGroupFreq", signature = "x",
    function(x) standardGeneric("tagSNVGroupFreq")
)

setGeneric("tagSNVGroupFreq<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVGroupFreq<-")
)



setGeneric("tagSNVChange", signature = "x",
    function(x) standardGeneric("tagSNVChange")
)

setGeneric("tagSNVChange<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVChange<-")
)


setGeneric("tagSNVsPerIndividual", signature = "x",
    function(x) standardGeneric("tagSNVsPerIndividual")
)

setGeneric("tagSNVsPerIndividual<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVsPerIndividual<-")
)


setGeneric("individualPerTagSNV", signature = "x",
    function(x) standardGeneric("individualPerTagSNV")
)

setGeneric("individualPerTagSNV<-", signature = c("x", "value"),
    function(x, value) standardGeneric("individualPerTagSNV<-")
)


setGeneric("tagSNVAnno", signature = "x",
    function(x) standardGeneric("tagSNVAnno")
)

setGeneric("tagSNVAnno<-", signature = c("x", "value"),
    function(x, value) standardGeneric("tagSNVAnno<-")
)




## class HaploClusterList


setGeneric("haploClusters", signature = "x",
    function(x) standardGeneric("haploClusters")
)

setGeneric("haploClusters<-", signature = c("x", "value"),
    function(x, value) standardGeneric("haploClusters<-")
)

setGeneric("lengthList", signature = "x",
    function(x) standardGeneric("lengthList")
)

setGeneric("lengthList<-", signature = c("x", "value"),
    function(x, value) standardGeneric("lengthList<-")
)
setGeneric("statistics", signature = "x",
    function(x) standardGeneric("statistics")
)

setGeneric("statistics<-", signature = c("x", "value"),
    function(x, value) standardGeneric("statistics<-")
)


setGeneric("setAnnotation", signature = c("haploClusterList","filename"),
    function(haploClusterList,filename) standardGeneric("setAnnotation")
)

setGeneric("setStatistics", signature = c("haploClusterList"),
    function(haploClusterList) standardGeneric("setStatistics")
)

setGeneric("compareHaploClusterLists", signature = c("haploClusterList1","haploClusterList2","simv","pTagSNVs","pIndivid","minTagSNVs","minIndivid"),
    function(haploClusterList1,haploClusterList2,simv,pTagSNVs,pIndivid,minTagSNVs,minIndivid) standardGeneric("compareHaploClusterLists")
)


setGeneric("mergeHaploClusterLists", signature = c("haploClusterList1","haploClusterList2","clustHaploClustList"),
    function(haploClusterList1,haploClusterList2,clustHaploClustList) standardGeneric("mergeHaploClusterLists")
)


setGeneric("haploClusterList2excel", signature = c("haploClusterList","filename"),
    function(haploClusterList,filename) standardGeneric("haploClusterList2excel")
)

