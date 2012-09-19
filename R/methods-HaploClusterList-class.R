### ------------------------------
### HaploClusterList class methods
### ------------------------------


##
## Constructor
##

HaploClusterList <- function(haploClusters=list(),lengthList=0,statistics=list()) {
    new("HaploClusterList", haploClusters=haploClusters,lengthList=lengthList,statistics=statistics)

}

##
## Getters and setters
##


setMethod("haploClusters", "HaploClusterList",
    function(x)
    {
      slot(x, "haploClusters")
    }
)


setReplaceMethod("haploClusters", c("HaploClusterList", "list"),
    function(x, value)
    {
       slot(x, "haploClusters") <- value
       x
    }
)


setMethod("lengthList", "HaploClusterList",
    function(x)
    {
      slot(x, "lengthList")
    }
)


setReplaceMethod("lengthList", c("HaploClusterList", "numeric"),
    function(x, value)
    {
       slot(x, "lengthList") <- value
       x
    }
)


setMethod("statistics", "HaploClusterList",
    function(x)
    {
      slot(x, "statistics")
    }
)


setReplaceMethod("statistics", c("HaploClusterList", "list"),
    function(x, value)
    {
       slot(x, "statistics") <- value
       x
    }
)




##
## Subsetting
##


setMethod("[[", c(x="HaploClusterList", i="numeric",j="missing"),
    function(x, i, ...)
      {
        haploClusters(x)[[i]]
      }
)


setReplaceMethod("[[", c(x="HaploClusterList",i="numeric",j="missing",value="HaploCluster"),
    function(x, i, ..., value)
    {
        haploClusters(x)[[i]] <- value
        if(i>lengthList(x)) {
            lengthList(x) <- i
        }
        x

    }
)

setMethod("[", c(x="HaploClusterList", i="numeric",j="missing"),
    function(x, i, ...)
      {
        new("HaploClusterList", haploClusters=haploClusters(x)[i],lengthList=length(i),statistics=list())
      }
)


setReplaceMethod("[", c(x="HaploClusterList",i="numeric",j="missing",value="HaploClusterList"),
    function(x, i, ..., value)
    {
        haploClusters(x)[i] <- value
        if(i>lengthList(x)) {
            lengthList(x) <- i
        }
        x

    }
)




##
## Summary
##



setMethod("summary", "HaploClusterList",
function(object, ...)
{
    cat("\nAn object of class",class(object))

    cat("\nNumber of haplotype clusters: ",lengthList(object))
    cat("\nStatistics:\n")
    print(statistics(object))

})



##
## Plot
##

setMethod("plot",signature(x="HaploClusterList", y="missing"),
function(x,filename, ...) {

    require(fabia)

    if (missing(x)) {
        stop("List of haplotype clusters 'x' is missing. Stopped.")
    }

    if (missing(filename)) {
        stop("File name 'filename' for loading genotypes is missing. Stopped.")
    }

devAskNewPage(ask = TRUE)

for (i in 1:lengthList(x)) {


individ <- individuals(x[[i]])
tagSNV <- tagSNVs(x[[i]])
tagSNV <- as.integer(sort.int(as.integer(unique(tagSNV))))
tagSNVPositions <- tagSNVPositions(x[[i]])
chrom <- chromosome(x[[i]])


labels_ALL <- labelIndividuals(x[[i]])

Lout <- readSamplesSpfabia(X=filename,samples=individ,lowerB=0,upperB=1000.0)


plotHaplotypeCluster(Lout=Lout,tagSNV=list(tagSNV),physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),chrom=chrom,count=i,labelsNA=labels_ALL)

}


devAskNewPage(ask = FALSE)

}
)


setMethod("setAnnotation",signature(haploClusterList="HaploClusterList",filename="character"),
function(haploClusterList,filename) {

nohaploClusters <- lengthList(haploClusterList)


if (nohaploClusters > 0) {

annotatedtagSNVs <- read.table(file=filename,header = FALSE, sep = " ", quote = "",as.is = TRUE)
# tagSNV position at first column !!!

anno_pos <- paste(annotatedtagSNVs[,2],"_",annotatedtagSNVs[,1],sep="")



for (haploClusterC in 1:nohaploClusters)  {


vt <- haploClusterList[[haploClusterC]]

chromosome <- chromosome(vt)
tagSNVPositions <- tagSNVPositions(vt)
vt_pos <- paste(chromosome(vt),"_",tagSNVPositions(vt),sep="")
annotMatchT <- match(vt_pos,anno_pos,nomatch=0)
tagSNVMatch <- which(annotMatchT>0)
annotMatch <- annotMatchT[which(annotMatchT>0)]
matchPositionhaploClusters <- tagSNVPositions[tagSNVMatch]
for (ianno in 1:numbertagSNVs(vt)) {
    tagSNVAnno(vt)[ianno] <- list("-")
}
for (ianno in 1:length(annotMatch)) {
    tagSNVAnno(vt)[tagSNVMatch] <- as.list(annotatedtagSNVs[annotMatch,c(-1,-2)])
}

}
}

haploClusterList

}
)


setMethod("setStatistics",signature(haploClusterList="HaploClusterList"),
function(haploClusterList) {



avhaploClusterPos <- c()
avhaploClusterLengthSNV <- c()
avhaploClusterLength <- c()
avnoIndivid <- c()
avnoTagSNVs <- c()

avnoFreq <- c()
avnoGroupFreq <- c()
avnotagSNVChange <- c()
avnotagSNVsPerIndividual <- c()
avnoindividualPerTagSNV <- c()


nohaploClusters <- lengthList(haploClusterList)




if (nohaploClusters > 0) {

for (haploClusterC in 1:nohaploClusters)  {


vt <- haploClusterList[[haploClusterC]]

avhaploClusterPos <- c(avhaploClusterPos,haploClusterPos(vt))
avhaploClusterLengthSNV <- c(avhaploClusterLengthSNV,haploClusterLength(vt))
avhaploClusterLength <- c(avhaploClusterLength,(max(tagSNVPositions(vt))- min(tagSNVPositions(vt))))
avnoIndivid <- c(avnoIndivid,numberIndividuals(vt))
avnoTagSNVs <- c(avnoTagSNVs,numbertagSNVs(vt))

avnoFreq <- c(avnoFreq,tagSNVFreq(vt))
avnoGroupFreq <- c(avnoGroupFreq,tagSNVGroupFreq(vt))
avnotagSNVChange <- c(avnotagSNVChange,tagSNVChange(vt))
avnotagSNVsPerIndividual <- c(avnotagSNVsPerIndividual,tagSNVsPerIndividual(vt))
avnoindividualPerTagSNV <- c(avnoindividualPerTagSNV,individualPerTagSNV(vt))

}
}


avhaploClusterPosS <- summary(avhaploClusterPos)
avhaploClusterLengthSNVS <- summary(avhaploClusterLengthSNV)
avhaploClusterLengthS <- summary(avhaploClusterLength)
avnoIndividS <- summary(avnoIndivid)
avnoTagSNVsS <- summary(avnoTagSNVs)

avnoFreqS <- summary(avnoFreq)
avnoGroupFreqS <- summary(avnoGroupFreq)
avnotagSNVChangeS <- summary(avnotagSNVChange)
avnotagSNVsPerIndividualS <- summary(avnotagSNVsPerIndividual)
avnoindividualPerTagSNVS <- summary(avnoindividualPerTagSNV)



statistics(haploClusterList) <- list(avhaploClusterPosS=avhaploClusterPosS,avhaploClusterLengthSNVS=avhaploClusterLengthSNVS,avhaploClusterLengthS=avhaploClusterLengthS,avnoIndividS=avnoIndividS,avnoTagSNVsS=avnoTagSNVsS,avnoFreqS=avnoFreqS,avnoGroupFreqS=avnoGroupFreqS,avnotagSNVChangeS=avnotagSNVChangeS,avnotagSNVsPerIndividualS=avnotagSNVsPerIndividualS,avnoindividualPerTagSNVS=avnoindividualPerTagSNVS)

haploClusterList

}
)



setMethod("compareHaploClusterLists",signature(haploClusterList1="HaploClusterList",haploClusterList2="ANY",simv="character",pTagSNVs="ANY",pIndivid="ANY",minTagSNVs="numeric",minIndivid="numeric"),
function(haploClusterList1,haploClusterList2=NULL,simv="minD",pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=6,minIndivid=2) {

    if (missing(haploClusterList1)) {
        stop("List of haplotype clusters 'haploClusterList1' is missing. Stopped.")
    }

    if (!is.null(haploClusterList2)) {
     tagSNVs1 <- sapply(haploClusters(haploClusterList1),function(x) {tagSNVs(x)} , simplify=FALSE)
     tagSNVs2 <- sapply(haploClusters(haploClusterList2),function(x) {tagSNVs(x)}  , simplify=FALSE)
     individ1 <- sapply(haploClusters(haploClusterList1),function(x) {individuals(x)}  , simplify=FALSE)
     individ2 <- sapply(haploClusters(haploClusterList2),function(x) {individuals(x)}  , simplify=FALSE)

     tagSNVs <- c(tagSNVs1,tagSNVs2)
     individ <- c(individ1,individ2)

  } else {

     tagSNVs <- sapply(haploClusters(haploClusterList1),function(x) {tagSNVs(x)} , simplify=FALSE )
     individ <- sapply(haploClusters(haploClusterList1),function(x) {individuals(x)}  , simplify=FALSE)
  }

  l <- length(tagSNVs)

  if (l>1) {

  tagSNVsSim <- sapply(tagSNVs,function(x) { sapply(tagSNVs, function(y) {sim(x,y,simv,minInter=minTagSNVs)})} )

  if (!is.null(pTagSNVs)) {
    tagSNVsSim <- tagSNVsSim^pTagSNVs
  }

  individSim <- sapply(individ,function(x) { sapply(individ, function(y) {sim(x,y,simv,minInter=minIndivid)})} )

  if (!is.null(pIndivid)) {
      individSim <- individSim^pIndivid
  }
  oneM <- matrix(1,nrow=l,ncol=l)

  dist <- as.dist(oneM - tagSNVsSim*individSim)


  clust <- hclust(dist)

  } else {

      clust=NULL
  }

  return(clust)

}
)


setMethod("mergeHaploClusterLists",signature(haploClusterList1="HaploClusterList",haploClusterList2="ANY",clustHaploClustList="vector"),
function(haploClusterList1,haploClusterList2=NULL,clustHaploClustList) {

    if (missing(haploClusterList1)) {
        stop("List of haplotype clusters 'haploClusterList1' is missing. Stopped.")
    }

    if (missing(clustHaploClustList)) {
        stop("Vector 'clustHaploClustList' of new cluster membership for each old cluster is missing. Stopped.")
    }


  l <- lengthList(haploClusterList1)
  if (!is.null(haploClusterList2)) {
      l2 <- lengthList(haploClusterList2)
  } else {
      l2 <- 0
  }

  tabClust <- table(clustHaploClustList)
  cl <- length(tabClust)

  haploClusterListmerge <-  new("HaploClusterList", haploClusters=list(),lengthList=0,statistics=list())

  for (i in 1:cl) {


    tr1 <- which(clustHaploClustList==i)
    ltr1 <- length(tr1)

    vv <- list()

    p2 <- tr1[1]


    if (p2<=l) {
        vv <- haploClusterList1[[p2]]
    } else  {
       if ((l2>0)&&(p2<=l+l2)) {
          p2 <- p2-l
          vv <- haploClusterList2[[p2]]
        }
    }

    if (ltr1>1) {
    for (j in 2:ltr1) {

        doit <- FALSE
        p2 <- tr1[j]
        if (p2<=l) {
            vt <- haploClusterList1[[p2]]
            doit <- TRUE
        } else {
            if ((l2>0)&&(p2<=(l+l2))) {
                p2 <- p2-l
                vt <- haploClusterList2[[p2]]
                doit <- TRUE
            }
        }

        if (doit) {

            ri1 <- individuals(vt)
            ri2 <- tagSNVs(vt)

            ma1 <- match(ri1,individuals(vv))
            ma2 <- match(ri2,tagSNVs(vv))
            na1 <- which(ma1>0)
            na2 <- which(ma2>0)

            if (length(na1)>0) {
                ma1 <- ma1[na1]
                tagSNVsPerIndividual(vv)[ma1] <- pmax(tagSNVsPerIndividual(vv)[ma1],tagSNVsPerIndividual(vt)[na1])
            }

            if (length(na2)>0) {
                ma2 <- ma2[na2]
                individualPerTagSNV(vv)[ma2] <- pmax(individualPerTagSNV(vv)[ma2],individualPerTagSNV(vt)[na2])
            }


            ab1 <- 1:numberIndividuals(vt)
            ab2 <- 1:numbertagSNVs(vt)

            a1 <- ab1[-na1]
            a2 <- ab2[-na2]

            if (length(a1)>0) {
                individuals(vv) <- c(individuals(vv),individuals(vt)[a1])
                populationIndividuals(vv) <- c(populationIndividuals(vv),populationIndividuals(vt)[a1])
                idIndividuals(vv) <- c(idIndividuals(vv),idIndividuals(vt)[a1])
                labelIndividuals(vv) <- c(labelIndividuals(vv),labelIndividuals(vt)[a1])
                platformIndividuals(vv) <- c(platformIndividuals(vv),platformIndividuals(vt)[a1])
                tagSNVsPerIndividual(vv) <- c(tagSNVsPerIndividual(vv),tagSNVsPerIndividual(vt)[a1])
            }

            if (length(a2)>0) {
                tagSNVs(vv) <- c(tagSNVs(vv),tagSNVs(vt)[a2])
                tagSNVPositions(vv) <- c(tagSNVPositions(vv),tagSNVPositions(vt)[a2])
                tagSNVAlleles(vv) <- c(tagSNVAlleles(vv),tagSNVAlleles(vt)[a2])
                tagSNVNames(vv) <- c(tagSNVNames(vv),tagSNVNames(vt)[a2])
                tagSNVFreq(vv) <- c(tagSNVFreq(vv),tagSNVFreq(vt)[a2])
                tagSNVGroupFreq(vv) <- c(tagSNVGroupFreq(vv),tagSNVGroupFreq(vt)[a2])
                tagSNVChange(vv) <- c(tagSNVChange(vv),tagSNVChange(vt)[a2])
                individualPerTagSNV(vv) <- c(individualPerTagSNV(vv),individualPerTagSNV(vt)[a2])
                tagSNVAnno(vv) <- c(tagSNVAnno(vv),tagSNVAnno(vt)[a2])
            }

            coreClusterIndividuals(vv) <- c(coreClusterIndividuals(vv),coreClusterIndividuals(vt))
        }


    }



    }


    so8 <- sort(individuals(vv),index.return=TRUE)

    individuals(vv) <- individuals(vv)[so8$ix]
    populationIndividuals(vv) <- populationIndividuals(vv)[so8$ix]
    idIndividuals(vv) <- idIndividuals(vv)[so8$ix]
    labelIndividuals(vv) <- labelIndividuals(vv)[so8$ix]
    platformIndividuals(vv) <- platformIndividuals(vv)[so8$ix]
    tagSNVsPerIndividual(vv) <- tagSNVsPerIndividual(vv)[so8$ix]



    so9 <- sort(tagSNVs(vv),index.return=TRUE)

    tagSNVs(vv) <- tagSNVs(vv)[so9$ix]
    tagSNVPositions(vv) <- tagSNVPositions(vv)[so9$ix]
    tagSNVAlleles(vv) <- tagSNVAlleles(vv)[so9$ix]
    tagSNVNames(vv) <- tagSNVNames(vv)[so9$ix]
    tagSNVFreq(vv) <- tagSNVFreq(vv)[so9$ix]
    tagSNVGroupFreq(vv) <- tagSNVGroupFreq(vv)[so9$ix]
    tagSNVChange(vv) <- tagSNVChange(vv)[so9$ix]
    individualPerTagSNV(vv) <- individualPerTagSNV(vv)[so9$ix]
    tagSNVAnno(vv) <- tagSNVAnno(vv)[so9$ix]


    haploClusterPos(vv) <- round(median(tagSNVPositions(vv)))
    haploClusterLength(vv) <- max(tagSNVs(vv))-min(tagSNVs(vv))
    numberIndividuals(vv) <- length(individuals(vv))
    numbertagSNVs(vv) <- length(tagSNVs(vv))
    coreClusterIndividuals(vv) <- unique(coreClusterIndividuals(vv))






   ID(vv) <- i
   bicluster_id(vv) <- ltr1



   haploClusterListmerge[[i]] <- vv



}

 return(haploClusterListmerge)

}
)




setMethod("haploClusterList2excel",signature(haploClusterList="HaploClusterList",filename="character"),
function(haploClusterList,filename) {

    if (missing(haploClusterList)) {
        stop("Object 'haploClusterList' of class HaploClusterList is missing. Stopped.")
    }
    if (missing(filename)) {
        stop("File name 'filename' of the EXCEL output file is missing. Stopped.")
    }

out <- c("ID , ","bicluster ID , ","chromosome , ","haploClusterPos , ","haploClusterLength , ","numberIndividuals , ","numbertagSNVs , ","IndividualNumber , ","TagSNVNumber , ","populationIndividuals , ","idIndividuals , ","labelIndividuals , ","platformIndividuals , ","coreClusterIndividuals , ","tagSNVPositions , ","tagSNVAlleles , ","tagSNVNames , ","tagSNVFreq , ","tagSNVGroupFreq , ","tagSNVChange , ","tagSNVsPerIndividual , ","individualPerTagSNV , ","tagSNVAnno")


write.table(t(out),file=filename,sep=" ",quote = FALSE,row.names = FALSE,col.names = FALSE)


lhaploCluster <- lengthList(haploClusterList)
if (lhaploCluster>0)  {
for (i in 1:lhaploCluster) {

    haplotClusti <- haploClusterList[[i]]
    numbertagSNVs <- numbertagSNVs(haplotClusti)

    out <- c(ID(haplotClusti))
    out <- c(out,",",bicluster_id(haplotClusti))
    out <- c(out,",",chromosome(haplotClusti))
    out <- c(out,",",haploClusterPos(haplotClusti))
    out <- c(out,",",haploClusterLength(haplotClusti))
    out <- c(out,",",numberIndividuals(haplotClusti))
    out <- c(out,",",numbertagSNVs)
    out <- c(out,",",individuals(haplotClusti))
    out <- c(out,",",tagSNVs(haplotClusti))
    out <- c(out,",",populationIndividuals(haplotClusti))
    out <- c(out,",",idIndividuals(haplotClusti))
    out <- c(out,",",labelIndividuals(haplotClusti))
    out <- c(out,",",platformIndividuals(haplotClusti))
    out <- c(out,",",coreClusterIndividuals(haplotClusti))
    out <- c(out,",",tagSNVPositions(haplotClusti))
    out <- c(out,",",tagSNVAlleles(haplotClusti))
    out <- c(out,",",tagSNVNames(haplotClusti))
    out <- c(out,",",tagSNVFreq(haplotClusti))
    out <- c(out,",",tagSNVGroupFreq(haplotClusti))
    out <- c(out,",",tagSNVChange(haplotClusti))
    out <- c(out,",",tagSNVsPerIndividual(haplotClusti))
    out <- c(out,",",individualPerTagSNV(haplotClusti))

    annoT <- unlist(tagSNVAnno(haplotClusti))
    lele <- length(annoT)
    if (length(table(annoT))>1) {
        segs <- lele /numbertagSNVs
        nam <- names(annoT)
        if (length(nam)<lele) {
            nam <- rep(as.character(1:segs),numbertagSNVs)
        }
        lau <- 1
        outT <- ""
        for (tuz in 1:numbertagSNVs)
        {
            if (tuz>1) {
                outT <- paste(outT," || ",sep="")
            }
            for (taz in 1:segs)
                outT <- paste(outT,nam[lau],"=",annoT[lau]," ",sep="")
        }

        out <- c(out,",",outT)
    } else {

        out <- c(out,","," ")

    }
    write.table(t(out),file=filename,sep=" ",quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)

}
}

}
)


