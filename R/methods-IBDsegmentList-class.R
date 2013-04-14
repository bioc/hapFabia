### ------------------------------
### IBDsegmentList class methods
### ------------------------------


##
## Constructor
##

IBDsegmentList <- function(IBDsegments=list(),lengthList=0,statistics=list()) {
    new("IBDsegmentList", IBDsegments=IBDsegments,lengthList=lengthList,statistics=statistics)

}

##
## Getters and setters
##


setMethod("IBDsegments", "IBDsegmentList",
    function(x)
    {
      slot(x, "IBDsegments")
    }
)


setReplaceMethod("IBDsegments", c("IBDsegmentList", "list"),
    function(x, value)
    {
       slot(x, "IBDsegments") <- value
       x
    }
)


setMethod("lengthList", "IBDsegmentList",
    function(x)
    {
      slot(x, "lengthList")
    }
)


setReplaceMethod("lengthList", c("IBDsegmentList", "numeric"),
    function(x, value)
    {
       slot(x, "lengthList") <- value
       x
    }
)


setMethod("statistics", "IBDsegmentList",
    function(x)
    {
      slot(x, "statistics")
    }
)


setReplaceMethod("statistics", c("IBDsegmentList", "list"),
    function(x, value)
    {
       slot(x, "statistics") <- value
       x
    }
)




##
## Subsetting
##


setMethod("[[", c(x="IBDsegmentList", i="numeric",j="missing"),
    function(x, i, ...)
      {
        IBDsegments(x)[[i]]
      }
)


setReplaceMethod("[[", c(x="IBDsegmentList",i="numeric",j="missing",value="IBDsegment"),
    function(x, i, ..., value)
    {
        IBDsegments(x)[[i]] <- value
        if(i>lengthList(x)) {
            lengthList(x) <- i
        }
        x

    }
)

setMethod("[", c(x="IBDsegmentList", i="numeric",j="missing"),
    function(x, i, ...)
      {
        new("IBDsegmentList", IBDsegments=IBDsegments(x)[i],lengthList=length(i),statistics=list())
      }
)


setReplaceMethod("[", c(x="IBDsegmentList",i="numeric",j="missing",value="IBDsegmentList"),
    function(x, i, ..., value)
    {
        IBDsegments(x)[i] <- value
        if(i>lengthList(x)) {
            lengthList(x) <- i
        }
        x

    }
)




##
## Summary
##



setMethod("summary", "IBDsegmentList",
function(object, ...)
{
    cat("\nAn object of class",class(object))

    cat("\nNumber of IBD segments: ",lengthList(object))
    cat("\nStatistics:\n")
    print(statistics(object))

})



##
## Plot
##

setMethod("plot",signature(x="IBDsegmentList", y="missing"),
function(x,filename, ...) {

    require(fabia)

    if (missing(x)) {
        stop("List of IBD segments 'x' is missing. Stopped.")
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


plotIBDsegment(Lout=Lout,tagSNV=list(tagSNV),physPos=tagSNVPositions,colRamp=12,val=c(0.0,2.0,1.0),chrom=chrom,count=i,labelsNA=labels_ALL)

}


devAskNewPage(ask = FALSE)

}
)


setMethod("setAnnotation",signature(IBDsegmentList="IBDsegmentList",filename="character"),
function(IBDsegmentList,filename) {

noIBDsegments <- lengthList(IBDsegmentList)


if (noIBDsegments > 0) {

annotatedtagSNVs <- read.table(file=filename,header = FALSE, sep = " ", quote = "",as.is = TRUE)
# tagSNV position at first column !!!

anno_pos <- paste(annotatedtagSNVs[,2],"_",annotatedtagSNVs[,1],sep="")



for (IBDsegmentC in 1:noIBDsegments)  {


vt <- IBDsegmentList[[IBDsegmentC]]

chromosome <- chromosome(vt)
tagSNVPositions <- tagSNVPositions(vt)
vt_pos <- paste(chromosome(vt),"_",tagSNVPositions(vt),sep="")
annotMatchT <- match(vt_pos,anno_pos,nomatch=0)
tagSNVMatch <- which(annotMatchT>0)
annotMatch <- annotMatchT[which(annotMatchT>0)]
matchPositionIBDsegments <- tagSNVPositions[tagSNVMatch]
for (ianno in 1:numbertagSNVs(vt)) {
    tagSNVAnno(vt)[ianno] <- list("-")
}
for (ianno in 1:length(annotMatch)) {
    tagSNVAnno(vt)[tagSNVMatch] <- as.list(annotatedtagSNVs[annotMatch,c(-1,-2)])
}

}
}

IBDsegmentList

}
)


setMethod("setStatistics",signature(IBDsegmentList="IBDsegmentList"),
function(IBDsegmentList) {



avIBDsegmentPos <- c()
avIBDsegmentLengthSNV <- c()
avIBDsegmentLength <- c()
avnoIndivid <- c()
avnoTagSNVs <- c()

avnoFreq <- c()
avnoGroupFreq <- c()
avnotagSNVChange <- c()
avnotagSNVsPerIndividual <- c()
avnoindividualPerTagSNV <- c()


noIBDsegments <- lengthList(IBDsegmentList)




if (noIBDsegments > 0) {

for (IBDsegmentC in 1:noIBDsegments)  {


vt <- IBDsegmentList[[IBDsegmentC]]

avIBDsegmentPos <- c(avIBDsegmentPos,IBDsegmentPos(vt))
avIBDsegmentLengthSNV <- c(avIBDsegmentLengthSNV,IBDsegmentLength(vt))
avIBDsegmentLength <- c(avIBDsegmentLength,(max(tagSNVPositions(vt))- min(tagSNVPositions(vt))))
avnoIndivid <- c(avnoIndivid,numberIndividuals(vt))
avnoTagSNVs <- c(avnoTagSNVs,numbertagSNVs(vt))

avnoFreq <- c(avnoFreq,tagSNVFreq(vt))
avnoGroupFreq <- c(avnoGroupFreq,tagSNVGroupFreq(vt))
avnotagSNVChange <- c(avnotagSNVChange,tagSNVChange(vt))
avnotagSNVsPerIndividual <- c(avnotagSNVsPerIndividual,tagSNVsPerIndividual(vt))
avnoindividualPerTagSNV <- c(avnoindividualPerTagSNV,individualPerTagSNV(vt))

}
}


avIBDsegmentPosS <- summary(avIBDsegmentPos)
avIBDsegmentLengthSNVS <- summary(avIBDsegmentLengthSNV)
avIBDsegmentLengthS <- summary(avIBDsegmentLength)
avnoIndividS <- summary(avnoIndivid)
avnoTagSNVsS <- summary(avnoTagSNVs)

avnoFreqS <- summary(avnoFreq)
avnoGroupFreqS <- summary(avnoGroupFreq)
avnotagSNVChangeS <- summary(avnotagSNVChange)
avnotagSNVsPerIndividualS <- summary(avnotagSNVsPerIndividual)
avnoindividualPerTagSNVS <- summary(avnoindividualPerTagSNV)



statistics(IBDsegmentList) <- list(avIBDsegmentPosS=avIBDsegmentPosS,avIBDsegmentLengthSNVS=avIBDsegmentLengthSNVS,avIBDsegmentLengthS=avIBDsegmentLengthS,avnoIndividS=avnoIndividS,avnoTagSNVsS=avnoTagSNVsS,avnoFreqS=avnoFreqS,avnoGroupFreqS=avnoGroupFreqS,avnotagSNVChangeS=avnotagSNVChangeS,avnotagSNVsPerIndividualS=avnotagSNVsPerIndividualS,avnoindividualPerTagSNVS=avnoindividualPerTagSNVS)

IBDsegmentList

}
)



setMethod("compareIBDsegmentLists",signature(IBDsegmentList1="IBDsegmentList",IBDsegmentList2="ANY",simv="character",pTagSNVs="ANY",pIndivid="ANY",minTagSNVs="numeric",minIndivid="numeric"),
function(IBDsegmentList1,IBDsegmentList2=NULL,simv="minD",pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=6,minIndivid=2) {

    if (missing(IBDsegmentList1)) {
        stop("List of IBD segments 'IBDsegmentList1' is missing. Stopped.")
    }

    if (!is.null(IBDsegmentList2)) {
     tagSNVs1 <- sapply(IBDsegments(IBDsegmentList1),function(x) {tagSNVs(x)} , simplify=FALSE)
     tagSNVs2 <- sapply(IBDsegments(IBDsegmentList2),function(x) {tagSNVs(x)}  , simplify=FALSE)
     individ1 <- sapply(IBDsegments(IBDsegmentList1),function(x) {individuals(x)}  , simplify=FALSE)
     individ2 <- sapply(IBDsegments(IBDsegmentList2),function(x) {individuals(x)}  , simplify=FALSE)

     tagSNVs <- c(tagSNVs1,tagSNVs2)
     individ <- c(individ1,individ2)

  } else {

     tagSNVs <- sapply(IBDsegments(IBDsegmentList1),function(x) {tagSNVs(x)} , simplify=FALSE )
     individ <- sapply(IBDsegments(IBDsegmentList1),function(x) {individuals(x)}  , simplify=FALSE)
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


setMethod("mergeIBDsegmentLists",signature(IBDsegmentList1="IBDsegmentList",IBDsegmentList2="ANY",clustIBDsegmentList="vector"),
function(IBDsegmentList1,IBDsegmentList2=NULL,clustIBDsegmentList) {

    if (missing(IBDsegmentList1)) {
        stop("List of IBD segments 'IBDsegmentList1' is missing. Stopped.")
    }

    if (missing(clustIBDsegmentList)) {
        stop("Vector 'clustIBDsegmentList' of new cluster membership for each old cluster is missing. Stopped.")
    }


  l <- lengthList(IBDsegmentList1)
  if (!is.null(IBDsegmentList2)) {
      l2 <- lengthList(IBDsegmentList2)
  } else {
      l2 <- 0
  }

  tabClust <- table(clustIBDsegmentList)
  cl <- length(tabClust)

  IBDsegmentListmerge <-  new("IBDsegmentList", IBDsegments=list(),lengthList=0,statistics=list())

  for (i in 1:cl) {


    tr1 <- which(clustIBDsegmentList==i)
    ltr1 <- length(tr1)

    vv <- list()

    p2 <- tr1[1]


    if (p2<=l) {
        vv <- IBDsegmentList1[[p2]]
    } else  {
       if ((l2>0)&&(p2<=l+l2)) {
          p2 <- p2-l
          vv <- IBDsegmentList2[[p2]]
        }
    }

    if (ltr1>1) {
    for (j in 2:ltr1) {

        doit <- FALSE
        p2 <- tr1[j]
        if (p2<=l) {
            vt <- IBDsegmentList1[[p2]]
            doit <- TRUE
        } else {
            if ((l2>0)&&(p2<=(l+l2))) {
                p2 <- p2-l
                vt <- IBDsegmentList2[[p2]]
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


    IBDsegmentPos(vv) <- round(median(tagSNVPositions(vv)))
    IBDsegmentLength(vv) <- max(tagSNVPositions(vv))-min(tagSNVPositions(vv))
    numberIndividuals(vv) <- length(individuals(vv))
    numbertagSNVs(vv) <- length(tagSNVs(vv))
    coreClusterIndividuals(vv) <- unique(coreClusterIndividuals(vv))






   ID(vv) <- i
   bicluster_id(vv) <- ltr1



   IBDsegmentListmerge[[i]] <- vv



}

 return(IBDsegmentListmerge)

}
)




setMethod("IBDsegmentList2excel",signature(IBDsegmentList="IBDsegmentList",filename="character"),
function(IBDsegmentList,filename) {

    if (missing(IBDsegmentList)) {
        stop("Object 'IBDsegmentList' of class IBDsegmentList is missing. Stopped.")
    }
    if (missing(filename)) {
        stop("File name 'filename' of the EXCEL output file is missing. Stopped.")
    }

out <- c("ID , ","bicluster ID , ","chromosome , ","IBDsegmentPos , ","IBDsegmentLength , ","numberIndividuals , ","numbertagSNVs , ","IndividualNumber , ","TagSNVNumber , ","populationIndividuals , ","idIndividuals , ","labelIndividuals , ","platformIndividuals , ","coreClusterIndividuals , ","tagSNVPositions , ","tagSNVAlleles , ","tagSNVNames , ","tagSNVFreq , ","tagSNVGroupFreq , ","tagSNVChange , ","tagSNVsPerIndividual , ","individualPerTagSNV , ","tagSNVAnno")


write.table(t(out),file=filename,sep=" ",quote = FALSE,row.names = FALSE,col.names = FALSE)


lIBDsegment <- lengthList(IBDsegmentList)
if (lIBDsegment>0)  {
for (i in 1:lIBDsegment) {

    IBDtsegmenti <- IBDsegmentList[[i]]
    numbertagSNVs <- numbertagSNVs(IBDtsegmenti)

    out <- c(ID(IBDtsegmenti))
    out <- c(out,",",bicluster_id(IBDtsegmenti))
    out <- c(out,",",chromosome(IBDtsegmenti))
    out <- c(out,",",IBDsegmentPos(IBDtsegmenti))
    out <- c(out,",",IBDsegmentLength(IBDtsegmenti))
    out <- c(out,",",numberIndividuals(IBDtsegmenti))
    out <- c(out,",",numbertagSNVs)
    out <- c(out,",",individuals(IBDtsegmenti))
    out <- c(out,",",tagSNVs(IBDtsegmenti))
    out <- c(out,",",populationIndividuals(IBDtsegmenti))
    out <- c(out,",",idIndividuals(IBDtsegmenti))
    out <- c(out,",",labelIndividuals(IBDtsegmenti))
    out <- c(out,",",platformIndividuals(IBDtsegmenti))
    out <- c(out,",",coreClusterIndividuals(IBDtsegmenti))
    out <- c(out,",",tagSNVPositions(IBDtsegmenti))
    out <- c(out,",",tagSNVAlleles(IBDtsegmenti))
    out <- c(out,",",tagSNVNames(IBDtsegmenti))
    out <- c(out,",",tagSNVFreq(IBDtsegmenti))
    out <- c(out,",",tagSNVGroupFreq(IBDtsegmenti))
    out <- c(out,",",tagSNVChange(IBDtsegmenti))
    out <- c(out,",",tagSNVsPerIndividual(IBDtsegmenti))
    out <- c(out,",",individualPerTagSNV(IBDtsegmenti))

    annoT <- unlist(tagSNVAnno(IBDtsegmenti))
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


