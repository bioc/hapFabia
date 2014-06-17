makePipelineFile  <- function(fileName,shiftSize=5000,intervalSize=10000,haplotypes=FALSE,dosage=FALSE) {
sourcefile <- file.path(system.file(package="hapFabia"),"pipeline/pipeline.R")

source1 <- readLines(sourcefile)

write(paste("#####define segments, overlap, filename #######",sep=""),file="pipeline.R",ncolumns=100)
write(paste("shiftSize <- ",shiftSize,sep=""),file="pipeline.R",append=TRUE,ncolumns=100)
write(paste("intervalSize <- ",intervalSize,sep=""),file="pipeline.R",append=TRUE,ncolumns=100)
write(paste("fileName <- \'",fileName,"\' # without type",sep=""),file="pipeline.R",append=TRUE,ncolumns=100)
if (haplotypes) {
    write("haplotypes <- TRUE",file="pipeline.R",append=TRUE,ncolumns=100)
} else {
    write("haplotypes <- FALSE",file="pipeline.R",append=TRUE,ncolumns=100)
}
if (dosage) {
    write("dosage <- TRUE",file="pipeline.R",append=TRUE,ncolumns=100)
} else {
    write("dosage <- FALSE",file="pipeline.R",append=TRUE,ncolumns=100)
}

write.table(source1,file="pipeline.R",append=TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)

}

hapFabiaVersion <- function() {
  version <- packageDescription("hapFabia",fields="Version")
  cat('\nhapFabia Package Version ', version, '\n')
  cat('\nCopyright (c) 2012 by Sepp Hochreiter\n\n')
}


findDenseRegions <- function(obs,p=0.9,inte=500,thres=11,off=0) {

    if (missing(obs)) {
        stop("Vector 'obs' for constructing histogram is missing. Stopped.")
    }

    off <- off%%inte
    bb <- seq(1-off,length(obs),inte)
    lb <- length(bb)
    bb <- c(bb,(bb[lb]+inte))
    a2 <- which(obs>quantile(obs,p))
    hh1 <- hist(a2,breaks=bb,plot = FALSE)

    a1 <- which(hh1$counts>thres)
    la <- length(a1)

    if (la>0) {

      b1 <- round(hh1$mid[a1])


      b2 <- rep(inte,la)

      b3 <- list()
      b4 <- rep(0,la)


      for (i in 1:la) {

        range <- (b1[i]-(inte%/%2)):(b1[i]+(inte%/%2))
        a3 <- intersect(range,a2)
        b3[[i]] <- sort(a3)
        b4[i] <- length(a3)
      }
      return(list(m=b1,l=b2,pos=b3,len=b4))
    } else {
      return(list(m=NULL,l=NULL,pos=NULL,len=NULL))
    }
}



matrixPlot <- function(x,range=NULL,yLabels=NULL,zlim=NULL,title=NULL,colRamp=12,grid=FALSE,pairs=FALSE,padj=NA,...){

     if (missing(x)) {
            stop("Data matrix x missing. Stopped.")
        }

        if (!is.matrix(x)) {
            x <- as.matrix(x)
        }

        if (!is.numeric(x)) {
            stop("Data matrix x must be numeric. Stopped.")
        }



     min <- min(x)
     max <- max(x)

      if( !is.null(zlim) ){
       min <- zlim[1]
       max <- zlim[2]
    }

    if( !is.null(yLabels) ){
       yLabels <- c(yLabels)
    } else {
       yLabels <- c(1:nrow(x))
   }

    if( !is.null(range) ){
      xLabels <- c(range)
    } else {
      xLabels <- c(1:ncol(x))
    }

     if( is.null(title) ){
     title <- c()
    }

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0

         seq00 <- seq(0,0,length=256)
         seq11 <- seq(1,1,length=256)

     if (colRamp<9) {
         seq01 <- seq(0,1,length=256)
         seq10 <- seq(1,0,length=256)
         seq051 <- seq(0.5,1,length=256)
         seq005 <- seq(0,0.5,length=256)

         dev1 <- cbind(seq01,seq01,seq10)
         dev2 <- cbind(seq01,seq051,seq10)
         dev3 <- cbind(seq01,seq005,seq10)
         blue2green <- cbind(seq00,seq01,seq10) # blue2green
         green2red <- cbind(seq01,seq10,seq00) # green2red
         blue2yellow <- cbind(seq01,seq01,seq10) # blue2yellow
         cyan2yellow <- cbind(seq01,seq11,seq10) # cyan2yellow
         magenta2green <- cbind(seq10,seq01,seq10) # magenta2green
     } else {
    if (colRamp==9) {



        seqP1 <- c(0,3,5,8,11,13,16,19,21,24,27,29,32,35,37,40,42,45,48,50,53,56,58,61,64,66,69,72,74,77,80,82,85,88,90,93,96,98,101,104,106,109,112,114,117,120,122,125,128,130,133,135,138,141,143,146,149,151,154,157,159,162,165,167,170,173,175,178,181,183,186,189,191,194,197,199,202,205,207,210,212,215,218,220,223,226,228,231,234,236,239,242,244,247,250,252) # 96 steps

        seqP2 <- rep(255,64)
        seqP3 <- seqP1[96:1]
        seqP4 <- rep(0,64)
        seqP5 <- seqP1[1:32]
        seqP6 <- seqP3[1:32]


        seqM1 <- c(seqP1,seqP2,seqP3)/255
        seqM2 <- c(seqP4,seqP1,seqP2,seqP6)/255
        seqM3 <- seqM1[256:1]
        matlab <- cbind(seqM1,seqM2,seqM3)
    } else {
        seq0055 <- seq(0,55/256,length=256)
        seqQ1 <- c(0,2,5,7,10,12,15,18,20,22,25,28,30,32,35,38,40,42,45,48,50,52,55,58,60,62,65,68,70,72,75,77,80,83,85,88,90,92,95,97,100,103,105,108,110,112,115,117,120,123,125,128,130,133,135,137,140,143,145,147,150,152,155,158,160,162,165,167,170,173,175,177,180,183,185,187,190,192,195,198,200,202,205,207,210,213,215,217,220,223,225,227,230,232,235,238,240,242,245,247,250,253) # 102 steps

        seqQ2 <- rep(255,52)
        seqQ3 <- seqQ1[102:77]
        seqQ4 <- rep(0,76)

        seqQ5 <- rep(0,26)
        seqQ6 <- (1:50)*5
        seqQ7 <- rep(255,104)
        seqQ8 <- seqQ6[50:1]

        seqQ9 <- seqQ3[26:1]
        seqQ10 <- rep(255,52)
        seqQ11 <- seqQ1[102:1]
        seqQ12 <- rep(0,76)

        seqL1 <- c(seqQ4,seqQ1,seqQ2,seqQ3)/255
        seqL2 <- c(seqQ5,seqQ6,seqQ7,seqQ8,seqQ5)/255
        seqL3 <- c(seqQ9,seqQ10,seqQ11,seqQ12)/255


        blue2green2red <- cbind(seqL1,seqL2,seqL3) # blue2green2red
        sepp1 <- cbind(seq11,seq11,seq11) - blue2green2red
        sepp2 <- sepp1
        sepp3 <- cbind(seq00,seq11,seq00) - blue2green2red
        sepp4 <- cbind(seq00,seq0055,seq00) - blue2green2red
        sepp5 <- cbind(seq11,seq00,seq11) - blue2green2red
        sepp6 <- sepp5

    }
}

    ColorRamp <- switch(colRamp,
    rgb(dev1),
    rgb(dev2),
    rgb(dev3),
    rgb(blue2green),
    rgb(green2red),
    rgb(blue2yellow),
    rgb(cyan2yellow),
    rgb(magenta2green),
    rgb(matlab),
    rgb(blue2green2red),
    rgb(blue2green2red),
    rgb(sepp1),
    rgb(sepp2),
    rgb(sepp3),
    rgb(sepp4),
    rgb(sepp5),
    rgb(sepp6),
      )

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
# par(mar = c(3,5,2.5,2))
# par(mar = c(3,5,2.5,1))
# Original: par(mar=c(5.1,4.1,4.1,2.1))
 marO <- par("mar")
 par(...)
 image(1:ncol(x), 1:nrow(x), t(x), col=ColorRamp, xlab="",
	ylab="", axes=FALSE, zlim=c(min,max))

 if( !is.null(title) ){
    title(main=title)
 }

xt <- c(1,axTicks(1))
lxt <- length(xt)

if (abs(xt[lxt]-ncol(x))>5) {
    xt <- c(xt,ncol(x))
}

if (grid) {
	grid(ncol(x), nrow(x), col="grey", lty = 1, lwd = 1)
}

if (pairs) {
	abline(h=seq(from =2.5, to=(nrow(x)-1), by=2), col="black")
}

axis(BELOW<-1, at=xt, labels=prettyNum(xLabels[xt], big.mark = ","), padj=padj)
axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1)

par(mar = marO)

layout(1)
}





plotIBDsegment <- function(Lout,tagSNV,physPos=NULL,colRamp=12,val=c(0.0,2.0,1.0),chrom="",count=0,labelsNA=NULL,prange=NULL,labelsNA1=c("model L"),grid=FALSE,pairs=FALSE,...) {

    if (missing(Lout)) {
        stop("Genotype data 'Lout' is missing. Stopped.")
    }

    if (missing(tagSNV)) {
        stop("tagSNVs 'tagSNV' are missing. Stopped.")
    }



n <- dim(Lout)[2]
n1 <- length(tagSNV)
doMis <- FALSE
if (n1>4) {
    n1 <- n1-1
    doMis <- TRUE
  }

if ((!is.null(prange))&&(!is.null(physPos))) {

if (length(prange)==2) {
ma <- length(tagSNV[[1]])

a1 <- which(physPos>=prange[1])
a2 <- which(physPos<=prange[2])

aa <- intersect(a1,a2)

if (length(aa)>1) {

  tagSNV <- tagSNV[[1]][aa]
  physPos <- physPos[aa]

}

}

}

if (is.null(physPos)) {
    physPos <- tagSNV*100
}


ma <- length(tagSNV[[1]])


range <- tagSNV[[1]][1]:tagSNV[[1]][ma]

lr <- tagSNV[[1]][ma]-tagSNV[[1]][1]+1

physRanget <- physPos[1]:physPos[ma]

lt <- length(physRanget)

ltr <- lt%/%(lr-1)

ltt <- c(1,(1:(lr-2))*ltr,lt)

physRange <- physRanget[ltt]

physLength <- round(lt/1000)

physPosM <- (physPos[ma] + physPos[1])%/%2


mat <- matrix(val[1],(n+n1+1),lr)
mat[(n+1),] <- val[2]

for (iB in 1:n1) {

  tagSNV_pos <- match(tagSNV[[iB]],range)

  mat[(n+1+iB),tagSNV_pos] <- val[3]
}
if (doMis) {
   tagSNV_pos <- match(tagSNV[[n1+1]],range)
   mat[(n+1+3),tagSNV_pos] <- 0.35*val[3]
}

tagSNV_pos <- match(tagSNV[[1]],range)

for (j in 1:n) {
    lla <- which(Lout[range,j]>=0.1)
    llb <- setdiff(1:lr,lla)
    llc <- intersect(tagSNV_pos,lla)
    lla <- setdiff(lla,llc)

    mat[j,lla] <- val[2]
    mat[j,llb] <- val[1]
    mat[j,llc] <- val[3]

}

if (is.null(labelsNA)) {
    labelsNA <- as.character(1:n)
}

if (length(labelsNA1)!=n1) {
  if (n1>1) {
  ade <- rep("",(n1-1))
  labelsNA1 = c("model L",ade)
} else {
  labelsNA1 = c("model L")
}
}

labelsNA <- c(labelsNA,"",labelsNA1)

if (count<=0) {
matrixPlot(mat,range=physRange,yLabels=labelsNA, title=paste("chr: ",chrom,"  ||  pos: ",prettyNum(physPosM, big.mark = ","),"  ||  length: ",physLength,"kbp  ||  #tagSNVs: ",ma,"  ||  #Individuals: ",n,sep=""),colRamp=colRamp, grid=grid, pairs=pairs, ...)
} else {
matrixPlot(mat,range=physRange,yLabels=labelsNA, title=paste("count: ",count,"  ||  chr: ",chrom,"  ||  pos: ",prettyNum(physPosM, big.mark = ","),"  ||  length: ",physLength,"kbp  ||  #tagSNVs: ",ma,"  ||  #Individuals: ",n,sep=""),colRamp=colRamp, grid=grid, pairs=pairs, ...)
}


}





sim <- function(x,y,simv="minD",minInter=2) {

    if ((missing(x))||(missing(y))) {
        stop("'x' or 'y' is missing. Stopped.")
    }

  intervv1 <- length(intersect(x,y))

  if (intervv1>=minInter) {
     erg <- switch(simv,
       jaccard = intervv1/length(union(x,y)),
       dice = 2*intervv1/(length(x)+length(y)),
       minD = intervv1/ min(length(x),length(y)),
       maxD = intervv1/ max(length(x),length(y))
     )
   } else {
     erg <- 0
   }

  return(erg)

}



iterateIntervals <- function(startRun=1,endRun,shift=5000,intervalSize=10000,annotationFile=NULL,fileName,prefixPath="",sparseMatrixPostfix="_mat",annotPostfix="_annot.txt",individualsPostfix="_individuals.txt",individuals=0,lowerBP=0,upperBP=0.05,p=10,iter=40,quant=0.01,eps=1e-5,alpha=0.03,cyc=50,non_negative=1,write_file=0,norm=0,lap=100.0,IBDsegmentLength=50,Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.03,haplotypes=FALSE,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100)
{


labelsA <- c()
annot <-  c()

ina <- as.numeric(readLines(paste(prefixPath,fileName,sparseMatrixPostfix,".txt",sep=""),n=2))
individualsN <- ina[1]
snvs <- ina[2]

save(individualsN,snvs,file=paste(fileName,"_All",".Rda",sep=""))


for (posAll in startRun:endRun)  {

start <- (posAll-1)*shift
end <- start + intervalSize

if (end > snvs)  {

  end <- snvs
}


pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")


if (is.null(annotationFile)) {
    labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
    if (length(labelsAA[,2])<2) {
        if (haplotypes) {
            lA <- as.vector(unlist(rbind(1:individualsN,1:individualsN)))
        } else {
            lA <- as.vector(1:individualsN)
        }
    } else {
        if (haplotypes) {
            lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
        } else {
            lA <- as.vector(labelsAA[,2])
        }
    }
    indiA <-  cbind(as.character(lA),as.character(lA),as.character(lA),as.character(lA))
} else {
    indit <- read.table(annotationFile, header = FALSE, sep = "\t", quote = "",as.is=TRUE)
    lind <- length(indit)
    if (lind<4) {
        inditA <- read.table(annotationFile, header = FALSE, sep = " ", quote = "",as.is=TRUE)
        if (length(inditA)>lind) {
            indit <- inditA
            lind <- length(inditA)
        }
    }

    if (haplotypes) {
        if (lind>0) {
            indi1 <- as.vector(unlist(rbind(indit[,1],indit[,1]))) # because haplotypes individuals are doubled
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(unlist(rbind(1:individualsN,1:individualsN)))
            } else {
                lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
            }
            indi1 <- as.character(lA)
        }
        if (lind>1) {
            indi2 <- as.vector(unlist(rbind(indit[,2],indit[,2])))
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(unlist(rbind(1:individualsN,1:individualsN)))
            } else {
                lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
            }
            indi2 <- as.character(lA)
        }
        if (lind>2) {
            indi3 <- as.vector(unlist(rbind(indit[,3],indit[,3])))
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(unlist(rbind(1:individualsN,1:individualsN)))
            } else {
                lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
            }
            indi3 <- as.character(lA)
        }
        if (lind>3) {
            indi4 <- as.vector(unlist(rbind(indit[,4],indit[,4])))
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(unlist(rbind(1:individualsN,1:individualsN)))
            } else {
                lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
            }
            indi4 <- as.character(lA)
        }
     } else {
        if (lind>0) {
            indi1 <- as.vector(indit[,1])
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(1:individualsN)
            } else {
                lA <- as.vector(labelsAA[,2])
            }
            indi1 <- as.character(lA)
        }
        if (lind>1) {
            indi2 <- as.vector(indit[,2])
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(1:individualsN)
            } else {
                lA <- as.vector(labelsAA[,2])
            }
            indi2 <- as.character(lA)
        }
        if (lind>2) {
            indi3 <- as.vector(indit[,3])
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(1:individualsN)
            } else {
                lA <- as.vector(labelsAA[,2])
            }
            indi3 <- as.character(lA)
        }
        if (lind>3) {
            indi4 <- as.vector(indit[,4])
        } else {
            labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
            if (length(labelsAA[,2])<2) {
                lA <- as.vector(1:individualsN)
            } else {
                lA <- as.vector(labelsAA[,2])
            }
            indi4 <- as.character(lA)
        }
    }
    indi1 <- gsub(",",";",indi1)
    indi2 <- gsub(",",";",indi2)
    indi3 <- gsub(",",";",indi3)
    indi4 <- gsub(",",";",indi4)
    indiA <- cbind(indi1,indi2,indi3,indi4)
    labelsA <- indiA
}

# Here other labels might be possible:
# labelsA[,i]: 1=id 2=subPopulation 3=population 4 = platform

#indit <- read.table("phase1_integrated_calls.20101123.ALL.panel", header = FALSE, sep = "\t", quote = "",as.is=TRUE)
#indi1 <- as.vector(unlist(rbind(indit[,1],indit[,1]))) # because haplotypes individuals are doubled
#indi2 <- as.vector(unlist(rbind(indit[,2],indit[,2])))
#indi3 <- as.vector(unlist(rbind(indit[,3],indit[,3])))
#indi4 <- as.vector(unlist(rbind(indit[,4],indit[,4])))
#indi4 <- gsub(",",";",indi4)
#indiA <- cbind(indi1,indi2,indi3,indi4)

resHapFabia <- hapFabia(fileName=fileName,prefixPath=prefixPath,sparseMatrixPostfix=sparseMatrixPostfix,annotPostfix=annotPostfix,individualsPostfix=individualsPostfix,labelsA=labelsA,pRange=pRange,individuals=individuals,lowerBP=lowerBP,upperBP=upperBP,p=p,iter=iter,quant=quant,eps=eps,alpha=alpha,cyc=cyc,non_negative=non_negative,write_file=write_file,norm=norm,lap=lap,IBDsegmentLength=IBDsegmentLength,Lt = Lt,Zt = Zt,thresCount=thresCount,mintagSNVsFactor=mintagSNVsFactor,pMAF=pMAF,haplotypes=haplotypes,cut=cut,procMinIndivids=procMinIndivids,thresPrune=thresPrune,simv=simv,minTagSNVs=minTagSNVs,minIndivid=minIndivid,avSNVsDist=avSNVsDist,SNVclusterLength=SNVclusterLength)


IBDsegmentList2excel(resHapFabia$mergedIBDsegmentList,paste(fileName,pRange,".csv",sep=""))

save(resHapFabia,annot,file=paste(fileName,pRange,"_resAnno",".Rda",sep=""))

}


}



identifyDuplicates <- function(fileName,startRun=1,endRun,shift=5000,intervalSize=10000) {


snvs <- 0
resHapFabia <- list()
labelsA <- c()
annot <-  c()

load(file=paste(fileName,"_All",".Rda",sep=""))

avIBDsegmentLength <- list()
avIBDsegmentPos <- list()
count <- 0

for (posAll in startRun:endRun)  {

start <- (posAll-1)*shift
end <- start + intervalSize

if (end > snvs)  {

  end <- snvs
}


pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")


load(file=paste(fileName,pRange,"_resAnno.Rda",sep=""))

mergedIBDsegmentList <- resHapFabia$mergedIBDsegmentList

noIBDsegments <- lengthList(mergedIBDsegmentList)


if (noIBDsegments>0)  {

count <- count + noIBDsegments

avIBDsegmentPos[[posAll]] <-  sapply(IBDsegments(mergedIBDsegmentList),function(x) {IBDsegmentPos(x)}  , simplify=FALSE)

avIBDsegmentLength[[posAll]] <-  sapply(IBDsegments(mergedIBDsegmentList),function(x) {IBDsegmentLength(x)}  , simplify=FALSE)

}
}


IBDsegmentPos <- unlist(avIBDsegmentPos)
IBDsegmentLength <- unlist(avIBDsegmentLength)

IBDsegmentSim <- cbind(IBDsegmentPos,IBDsegmentLength)

dups <- duplicated(IBDsegmentSim)

un <- which(dups==FALSE)




#### enumerate counts ####




allCount <- 0
allCount1 <- 0

resD <- list()
resDA <- list()


for (posAll in startRun:endRun)  {

start <- (posAll-1)*shift
end <- start + intervalSize

if (end > snvs)  {

  end <- snvs
}



pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")

load(file=paste(fileName,pRange,"_resAnno.Rda",sep=""))

mergedIBDsegmentList <- resHapFabia$mergedIBDsegmentList

noIBDsegments <- lengthList(mergedIBDsegmentList)

if (noIBDsegments>0)  {

for (IBDsegmentC in 1:noIBDsegments)  {

allCount <- allCount + 1

resD1 <- c(allCount,IBDsegmentC,posAll)
resDA[[allCount]] <- resD1

if (!dups[allCount]) {

allCount1 <- allCount1 + 1

resD1 <- c(allCount1,allCount,IBDsegmentC,posAll)
resD[[allCount1]] <- resD1



}


}

}


}


if ( length(resD) >0 ) {
rr <- unlist(resD)
l <-4

countsA1 <- matrix(rr,nrow=allCount1,ncol=l,byrow=TRUE)

colnames(countsA1) <- c("allCount1","allCount","IBDsegmentC","posAll")

rr <- unlist(resDA)
l <-3

countsA2 <- matrix(rr,nrow=allCount,ncol=l,byrow=TRUE)

colnames(countsA2) <- c("allCount","IBDsegmentC","posAll")

} else {

dups <- FALSE
un <- 0
countsA1 <- c(0,0,0,0)
dim(countsA1) <- c(1,4)
colnames(countsA1) <- c("allCount1","allCount","IBDsegmentC","posAll")
countsA2 <- c(0,0,0)
dim(countsA2) <- c(1,3)
colnames(countsA2) <- c("allCount","IBDsegmentC","posAll")

}

# write.table(countsA1,file=paste("countsA1",runIndex,".txt",sep=""))
# write.table(countsA2,file=paste("countsA2",runIndex,".txt",sep=""))



save(dups,un,countsA1,countsA2,file=paste("dups.Rda",sep=""))


}


analyzeIBDsegments <- function(fileName,runIndex="",annotPostfix="_annot.txt",startRun=1,endRun,shift=5000,intervalSize=10000) {

resHapFabia <- list()
countsA2 <- c()
snvs <- 0
mergedIBDsegmentList <- list()
dups <- c()
load(file=paste(fileName,"_All",".Rda",sep=""))
load(file=paste("dups.Rda",sep=""))

if ((startRun>1)&&(length(countsA2[,3])>1)) {
  tzz <- countsA2[which(countsA2[,3]<startRun),]
  offC <- max(tzz[,1])
} else {
  offC <-  0
}


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



allCount <- 0
allCount1 <- 0

for (posAll in startRun:endRun)  {

start <- (posAll-1)*shift
end <- start + intervalSize

if (end > snvs)  {

  end <- snvs
}



pRange <- paste("_",format(start,scientific=FALSE),"_",format(end,scientific=FALSE),sep="")

load(file=paste(fileName,pRange,"_resAnno.Rda",sep=""))

mergedIBDsegmentList <- resHapFabia$mergedIBDsegmentList

noIBDsegments <- lengthList(mergedIBDsegmentList)




if (noIBDsegments > 0) {

for (IBDsegmentC in 1:noIBDsegments)  {

allCount <- allCount + 1

if (!dups[allCount+offC]) {

allCount1 <- allCount1 + 1

vt <- mergedIBDsegmentList[[IBDsegmentC]]

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



}

}

noIBDsegments <- allCount1


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



save(startRun,endRun,noIBDsegments,avIBDsegmentPos,avIBDsegmentLengthSNV,avIBDsegmentLength,avnoIndivid,avnoTagSNVs,avnoFreq,avnoGroupFreq,avnotagSNVChange,avnotagSNVsPerIndividual,avnoindividualPerTagSNV,avIBDsegmentPosS,avIBDsegmentLengthSNVS,avIBDsegmentLengthS,avnoIndividS,avnoTagSNVsS,avnoFreqS,avnoGroupFreqS,avnotagSNVChangeS,avnotagSNVsPerIndividualS,avnoindividualPerTagSNVS,file=paste("analyzeResult",runIndex,".Rda",sep=""))

return(list(startRun=startRun,endRun=endRun,noIBDsegments=noIBDsegments,avIBDsegmentPos=avIBDsegmentPos,avIBDsegmentLengthSNV=avIBDsegmentLengthSNV,avIBDsegmentLength=avIBDsegmentLength,avnoIndivid=avnoIndivid,avnoTagSNVs=avnoTagSNVs,avnoFreq=avnoFreq,avnoGroupFreq=avnoGroupFreq,avnotagSNVChange=avnotagSNVChange,avnotagSNVsPerIndividual=avnotagSNVsPerIndividual,avnoindividualPerTagSNV=avnoindividualPerTagSNV,avIBDsegmentPosS=avIBDsegmentPosS,avIBDsegmentLengthSNVS=avIBDsegmentLengthSNVS,avIBDsegmentLengthS=avIBDsegmentLengthS,avnoIndividS=avnoIndividS,avnoTagSNVsS=avnoTagSNVsS,avnoFreqS=avnoFreqS,avnoGroupFreqS=avnoGroupFreqS,avnotagSNVChangeS=avnotagSNVChangeS,avnotagSNVsPerIndividualS=avnotagSNVsPerIndividualS,avnoindividualPerTagSNVS=avnoindividualPerTagSNVS))



}


hapFabia <- function(fileName,prefixPath="",sparseMatrixPostfix="_mat",annotPostfix="_annot.txt",individualsPostfix="_individuals.txt",labelsA=NULL,pRange="",individuals=0,lowerBP=0,upperBP=0.05,p=10,iter=40,quant=0.01,eps=1e-5,alpha=0.03,cyc=50,non_negative=1,write_file=0,norm=0,lap=100.0,IBDsegmentLength=50,Lt = 0.1,Zt = 0.2,thresCount=1e-5,mintagSNVsFactor=3/4,pMAF=0.03,haplotypes=FALSE,cut=0.8,procMinIndivids=0.1,thresPrune=1e-3,simv="minD",minTagSNVs=6,minIndivid=2,avSNVsDist=100,SNVclusterLength=100) {

# fileName:            the file name of the sparse matrix in sparse format.
# prefixPath:          path of the data file
# sparseMatrixPostfix: postfix string for the sparse matrix
# annotPostfix:        postfix string for the annotation file
# labelsA:             individual names as matrix individuals x 4
# prange:              for intervals indicates its range
# individuals:             vector of individuals which should be included into the analysis; default = 0 (all individuals)
# lowerBP:             lower bound for filtering the inputs columns, minimal MAF (however more than one occurence to remove private SNVs)
# upperBP:             Upper bound for filtering the inputs columns, minimal MAF
# p:                   no biclusters per iteration
# alpha:               sparseness loadings; default = 0.03
# iter:                number iterations
# quant:               percentage of Ls to remove in each iteration
# eps:                 lower bound for variational parameter lapla; default: 1e-5
# cyc:                 number of iterations; default = 50
# non_negative:        Non-negative factors and loadings if non_negative; default = 1 (yes).
# write_file:          results are written to files (L in sparse format), default = 0 (not written).
# norm:                data normalization; default = 1 (no normalization).
# lap:                 minimal value of the variational parameter; default = 100.0.
# IBDsegmentLength:           IBD segment length in kbp
# Lt:                  percentage of largest Ls to consider for IBD segment extraction
# Zt:                  percentage of largest Zs to consider for IBD segment extraction
# thresCount:          p-value of random histogram hit, default 1e-5
# mintagSNVsFactor:       percentage of segments overlap in IBD segments; 1/2 for large to 3/4 for for small intervals
# pMAF:                averaged and corrected minor allele frequency
# haplotypes:          haplotypes = phased genotypes -> two chromosomes per individual

message("                      ")
message("                      ")
message("Running hapFabia with:")
message("   Prefix string for file name of data files --------- : ",fileName)
message("   Path of data files ---------------------------------: ",prefixPath)
if (haplotypes) {
message("   Data consists of phased genotypes (haplotypes) -----")
} else {
message("   Data consists of unphased genotypes ----------------")
}
message("   Postfix string for file in sparse matrix format ----: ",sparseMatrixPostfix)
message("   Postfix string for file containing individual names : ",individualsPostfix)
if (is.null(labelsA)) {
message("   Individuals annotation is not supplied -------------")
} else {
message("   Individuals annotation is supplied -----------------")
}
message("   String indicating the interval that is analyzed ----: ",pRange)
if (length(individuals)>1) {
message("   Number of individuals included into the analysis ---: ",length(individuals))
} else {
message("   All individuals are included into the analysis -----: 0 = all individuals")
}
message("   Lower bound on MAF but more than one occurence -----: ",lowerBP  )
message("   Upper bound on MAF ---------------------------------: ",upperBP)
message("   Number of biclusters per iteration -----------------: ",p)
message("   Sparseness coefficient of the loadings -------------: ",alpha)
message("   Number of iterations -------------------------------: ",iter)
message("   Percentage of Ls to remove after each iteration ----: ",quant)
message("   Lower bound for variational parameter lapla --------: ",eps)
message("   Number of cycles per iteration ---------------------: ",cyc)
message("   Non-negative factors and loadings ------------------: ",non_negative)
message("   Results written to files ---------------------------: ",write_file)
message("   Data normalized ------------------------------------: ",norm)
message("   Minimal value of the variational parameter ---------: ",lap)
if (IBDsegmentLength>0) {
message("   IBD segment length in kbp --------------------------: ",IBDsegmentLength)
message("   Average distance between SNVs in bases -------------: ",avSNVsDist)
} else {
if (SNVclusterLength>0)
{
message("   Number of SNVs in histogram bin --------------------: ",SNVclusterLength)
} else {
message("   Number of SNVs in histogram bin --------------------: 100")
}
}
message("   % largest Ls for IBD segment extraction ------------: ",Lt)
message("   % largest Zs for IBD segment extraction ------------: ",Zt)
message("   p-value threshold for random histogram counts ------: ",thresCount)
message("   Min. % of segments overlap in IBD segments ---------: ",mintagSNVsFactor)
message("   Averaged and corrected minor allele frequency ------: ",pMAF)
message("   Cutoff for merging IBD segments --------------------: ",cut)
message("   % of cluster individuals a tagSNV must tag ---------: ",procMinIndivids)
message("   Threshold for pruning border tagSNVs ---------------: ",thresPrune)
message("   Similarity measure for merging clusters ------------: ",simv)
message("   Minimum matching tagSNVs for cluster similarity ----: ",minTagSNVs)
message("   Minimum matching individuals for cluster similarity : ",minIndivid)
message("                      ")
message("                      ")



# Compute internal parameters

# ps: quantile above which to consider Ls
ps <- 1 - Lt
# psZ: quantile above which to consider Zs
psZ <- 1 - Zt


# compute histogram length
if (IBDsegmentLength>0) {
    if (avSNVsDist>0) {
        inteA <- IBDsegmentLength*(1000/avSNVsDist)
    } else {
        inteA <- IBDsegmentLength*10
    }
} else {
    if (SNVclusterLength>0)
    {
       inteA <- SNVclusterLength
    } else {
       inteA <- 100
    }
}

ina <- as.numeric(readLines(paste(prefixPath,fileName,pRange,sparseMatrixPostfix,".txt",sep=""),n=2))
if (length(individuals)>1) {
    individualsN <- length(individuals)
} else {
    individualsN <- ina[1]
    individuals <- 0
}
snvs <- ina[2]

upperBindivid=upperBP*individualsN # remove common SNVs
lowerBindivid=max(1.5,lowerBP*individualsN) # remove private SNVs

kk <- 1
while ((snvs/inteA)*choose(individualsN,2)*(1-pbinom(kk,inteA,pMAF*pMAF))>thresCount) { kk <- kk +1}

thresA <- kk
mintagSNVs <- round(mintagSNVsFactor*thresA)

# End Compute internal parameters



# Fabia call
res <- spfabia(X=paste(prefixPath,fileName,pRange,sparseMatrixPostfix,sep=""),p=p,alpha=alpha,cyc=cyc,non_negative=non_negative,write_file=write_file,norm=norm,lap=lap,samples=individuals,iter=iter,quant=quant,lowerB=lowerBindivid,upperB=upperBindivid,eps=eps)


# Load individuals to Ls of interest: load minor alleles of the Ls

sPF <- samplesPerFeature(X=paste(prefixPath,fileName,pRange,sparseMatrixPostfix,sep=""),samples=individuals,lowerB=lowerBindivid,upperB=upperBindivid)

if (nchar(annotPostfix)>0) {

#annot[[1]] <- chromosome
#annot[[2]] <- phys. position
#annot[[3]] <- snvNames
#annot[[4]] <- snvMajor
#annot[[5]] <- snvMinor
#annot[[6]] <- quality
#annot[[7]] <- pass
#annot[[8]] <- info of vcf file
#annot[[9]] <- fields in vcf file
#annot[[10]] <- frequency
#annot[[11]] <- 1 = changed if major allele is actually minor allele otherwise 0


    annot <- read.table(paste(prefixPath,fileName,pRange,annotPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE,skip=2)


    for (i in 1:length(annot)) {

        annot[[i]] <- gsub(",",";",annot[[i]])

    }
    for (i in 1:length(annot)) {

        annot[[i]] <- gsub("TRUE","T",annot[[i]])

    }

    annot[[2]] <- as.numeric(annot[[2]]) # physical position
    annot[[10]] <- as.numeric(annot[[10]]) # SNV frequency
    annot[[11]] <- as.numeric(annot[[11]]) # changed

}

if (is.null(labelsA)) {
    labelsAA <- read.table(paste(prefixPath,fileName,individualsPostfix,sep=""),header = FALSE, sep = " ", quote = "",as.is = TRUE)
    if (haplotypes) {
        lA <- as.vector(unlist(rbind(labelsAA[,2],labelsAA[,2])))
    } else {
        lA <- as.vector(labelsAA[,2])
    }
    indiA <-  cbind(as.character(lA),as.character(lA),as.character(lA),as.character(lA))
} else {
    indiA <- labelsA
}

# Here other labels might be possible:
# labelsA[,i]: 1=id 2=subPopulation 3=population 4 = platform

#indit <- read.table("phase1_integrated_calls.20101123.ALL.panel", header = FALSE, sep = "\t", quote = "",as.is=TRUE)
#indi1 <- as.vector(unlist(rbind(indit[,1],indit[,1]))) # because haplotypes individuals are doubled
#indi2 <- as.vector(unlist(rbind(indit[,2],indit[,2])))
#indi3 <- as.vector(unlist(rbind(indit[,3],indit[,3])))
#indi4 <- as.vector(unlist(rbind(indit[,4],indit[,4])))
#indi4 <- gsub(",",";",indi4)
#indiA <- cbind(indi1,indi2,indi3,indi4)



# save fabia result
# save(res,sPF,annot,file=paste(fileName,pRange,"_res.Rda",sep=""))


# first haplotype extraction with offset 0

off1 <- 0


IBDsegmentList1 <- extractIBDsegments(res=res,sPF=sPF,annot=annot,chrom="",labelsA=indiA,ps=ps,psZ=psZ,inteA=inteA,thresA=thresA,mintagSNVs=mintagSNVs,off=off1,procMinIndivids=procMinIndivids,thresPrune=thresPrune)

# merge IBD segment lists


if ( lengthList(IBDsegmentList1) > 1) {
  comp <- compareIBDsegmentLists(IBDsegmentList1=IBDsegmentList1,IBDsegmentList2=NULL,simv=simv,pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=minTagSNVs,minIndivid=minIndivid)

  if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp,h=cut)
      mergedIBDsegmentList1 <- mergeIBDsegmentLists(IBDsegmentList1=IBDsegmentList1,IBDsegmentList2=NULL,clustIBDsegmentList=clustIBDsegmentList)
  } else {
      mergedIBDsegmentList1 <- IBDsegmentList1
  }
} else {

 mergedIBDsegmentList1 <- IBDsegmentList1
}


# second IBD extraction with offset half of the interval length


off2=inteA%/%2

IBDsegmentList2 <- extractIBDsegments(res=res,sPF=sPF,annot=annot,chrom="",labelsA=indiA,ps=ps,psZ=psZ,inteA=inteA,thresA=thresA,mintagSNVs=mintagSNVs,off=off2,procMinIndivids=procMinIndivids,thresPrune=thresPrune)

# merge IBD segments

if ( lengthList(IBDsegmentList2) > 1) {
  comp <- compareIBDsegmentLists(IBDsegmentList1=IBDsegmentList2,IBDsegmentList2=NULL,simv=simv,pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=minTagSNVs,minIndivid=minIndivid)

  if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp,h=cut)
      mergedIBDsegmentList2 <- mergeIBDsegmentLists(IBDsegmentList1=IBDsegmentList2,IBDsegmentList2=NULL,clustIBDsegmentList=clustIBDsegmentList)
  } else {
      mergedIBDsegmentList2 <- IBDsegmentList2
  }
} else {
  mergedIBDsegmentList2 <- IBDsegmentList2
}

# merge IBD segments of both extractions


if ( lengthList(mergedIBDsegmentList1) > 0) {

  if ( lengthList(mergedIBDsegmentList2) > 0) {

    comp12 <- compareIBDsegmentLists(IBDsegmentList1=mergedIBDsegmentList1,IBDsegmentList2=mergedIBDsegmentList2,simv=simv,pTagSNVs=NULL,pIndivid=NULL,minTagSNVs=minTagSNVs,minIndivid=minIndivid)

    if (!is.null(comp12)) {
        clustIBDsegmentList <- cutree(comp12,h=cut)
        mergedIBDsegmentList <- mergeIBDsegmentLists(IBDsegmentList1=mergedIBDsegmentList1,IBDsegmentList2=mergedIBDsegmentList2,clustIBDsegmentList=clustIBDsegmentList)
   } else {
        mergedIBDsegmentList <- mergedIBDsegmentList1
    }

} else {

mergedIBDsegmentList <- mergedIBDsegmentList1
}

} else {

  if ( lengthList(mergedIBDsegmentList2) > 0) {
    mergedIBDsegmentList <- mergedIBDsegmentList2
  } else {
    mergedIBDsegmentList <- mergedIBDsegmentList1
  }


}

mergedIBDsegmentList1 <- setStatistics(mergedIBDsegmentList1)
mergedIBDsegmentList2 <- setStatistics(mergedIBDsegmentList2)
mergedIBDsegmentList <- setStatistics(mergedIBDsegmentList)

# save(mergedIBDsegmentList,res,sPF,annot,IBDsegmentList1,IBDsegmentList2,mergedIBDsegmentList1,mergedIBDsegmentList2,file=paste(fileName,pRange,"_IBDsegmentList.Rda",sep=""))

return(list(mergedIBDsegmentList=mergedIBDsegmentList,res=res,sPF=sPF,annot=annot,IBDsegmentList1=IBDsegmentList1,IBDsegmentList2=IBDsegmentList2,mergedIBDsegmentList1=mergedIBDsegmentList1,mergedIBDsegmentList2=mergedIBDsegmentList2))

}





simulateIBDsegments <- function(fileprefix="dataSim",minruns=1,maxruns=100,snvs=10000,individualsN=100,avDistSnvs=100,avDistMinor=25,noImplanted=1,implanted=10,length=100,minors=20,mismatches=0,mismatchImplanted=0.5,overlap=50,noOverwrite=FALSE) {

## minruns: min number or runs
## maxruns: max number or runs
## snvs: number of SNVs
## individualsN: number of individuals
## avDistSnvs: average distance in bases between tagSNVs
## avDistMinor: average distance between minor alleles, thus 1/avDistMinor is the average MAF
## noImplanted: how many IBD segments
## implanted: in how many individuals is one IBD segment implanted
## length: length of IBD segments
## minors: number of tagSNVs
## mismatches: number of mismatches in hapltopye
## mismatchImplanted: mismatches in how many percent of IBD segments
## overlap: minimal IBD segment overlap between chromosomes
## noOverwrite: noOverwrite=TRUE ensures that an IBD segment is not superimposed by another IBD segment

write(paste("maxruns: ",maxruns,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),ncolumns=100)
write(paste("snvs: ",snvs,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("individualsN: ",individualsN,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("avDistSnvs: ",avDistSnvs,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("avDistMinor: ",avDistMinor,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("noImplanted: ",noImplanted,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("implanted: ",implanted,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("length: ",length,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("minors: ",minors,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
write(paste("mismatches: ",mismatches,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
if (mismatches>0) {
  write(paste("mismatchImplanted: ",mismatchImplanted,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)
}
write(paste("overlap: ",overlap,sep=""),file=paste(fileprefix,"Parameters.txt",sep=""),append=TRUE,ncolumns=100)


rateSnvs <- 1/avDistSnvs

rateMinor <- 1/avDistMinor



haploN <- 2*individualsN

dataA <- 1:individualsN


nucleotide <- c("A","T","C","G")




physPos <- vector("integer",snvs)
snvMajor <- vector("character",snvs)
snvMinor <- vector("character",snvs)

g2 <- 2*(1:snvs)
g1 <- g2 -1


s2 <- 2*dataA
s1 <- s2-1


id <- as.character(dataA)
zeroid <-  as.character(rep(0,individualsN))

famID <- id
dim(famID) <- c(individualsN,1)
ID <- id
dim(ID) <- c(individualsN,1)
patID <- zeroid
dim(patID) <- c(individualsN,1)
matID <- zeroid
dim(matID) <- c(individualsN,1)
sex <- zeroid
dim(sex) <- c(individualsN,1)
phen <- zeroid
dim(phen) <- c(individualsN,1)



snvNames <- as.character(1:snvs)

chr <- rep(1,snvs)

alleleI <- matrix(0,haploN,snvs)
alleleN <- matrix("N",haploN,snvs)


#################

for (run in minruns:maxruns) {

  write(noImplanted,file=paste(fileprefix,run,"Impl.txt",sep=""),ncolumns=100)


  distsSnvs <- rexp(snvs,rateSnvs)

  aa <- as.integer(0)
  for (i in 1:snvs) {
    aa <- aa + 1 + as.integer(round(distsSnvs[i]))
    physPos[i] <- aa

    mami <- sample(4,2)
    snvMajor[i] <- nucleotide[mami[1]]
    snvMinor[i] <- nucleotide[mami[2]]
  }


  for (s in 1:haploN) {
    for (i in 1:snvs) {
      if (runif(1)<rateMinor) {
        alleleI[s,i] <- 1
        alleleN[s,i] <- snvMinor[i]
      } else {
        alleleI[s,i] <- 0
        alleleN[s,i] <- snvMajor[i]
      }
    }
  }



  freq <- colSums(alleleI)/haploN
  pos <- physPos
  pos1 <- pos/1000000




  plinkCols <- cbind(famID,ID,patID,matID,sex,phen)

  plinkmap <- cbind(chr,snvNames,pos1,pos)

  plinkmap[,1] <- as.integer(plinkmap[,1])
  plinkmap[,2] <- as.integer(plinkmap[,2])
  plinkmap[,3] <- format(as.double(plinkmap[,3]),nsmall=6)
  plinkmap[,4] <- as.integer(plinkmap[,4])

  plinkLine <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE",snvNames)


  mcmc <- rbind(pos1,freq)

  initM <- matrix(0,individualsN,2*snvs)


  alleleIimp <- alleleI
  alleleNimp <- alleleN

  allinter <- list()
  allimpl <- list()

  for (noIm in 1:noImplanted) {

  notfoundIBDsegments <- TRUE
  while (notfoundIBDsegments) {


  start <- sample((snvs-length-1),1)

  end <- start+length-1

  inter <- start:end

  allinter[[noIm]] <- inter

  posMinor <- sample(length,minors)

  implantI <- rep(0,length)

  implantI[posMinor] <- 1

  implantN <- rep("N",length)

  ll <- 1
  for (i in inter) {
    if (implantI[ll]==0) {
      implantN[ll] <- snvMajor[i]
    } else {
      implantN[ll] <- snvMinor[i]
    }
    ll <- ll + 1
  }


  impl <- sample(haploN,implanted)

  allimpl[[noIm]] <- impl

  notfoundIBDsegments <- FALSE
  if ((noIm>1)&&(noOverwrite)) {
    for (no2Im in 1:(noIm-1)) {
      if ((length(intersect(impl,allimpl[[no2Im]]))>0)&&(length(intersect(inter,allinter[[no2Im]]))>0)) {
        notfoundIBDsegments <- TRUE
      }
    }
  }


  }


  impl=sort(impl)

  implG <- (impl+1)%/%2
  write(impl,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  write(implG,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  write(inter,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  write(implantI,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  write(implantN,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)

  intervalImp <- matrix(0,implanted,4)


  for (i in 1:implanted) {


    startI <- 1+sample((length-overlap)%/%2,1)
    endI <- length - 1 - sample((length-overlap)%/%2,1)
    interI <- inter[startI:endI]
    alleleIimp[impl[i],interI] <- implantI[startI:endI]
    alleleNimp[impl[i],interI] <- implantN[startI:endI]

    interImp <- c(inter[startI],inter[endI],startI,endI)
    intervalImp[i,] <-  interImp

    write(interImp,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  }

  if (mismatches > 1) {

    posMismatcht <- 1+ sample((length-2),mismatches)

    posMismatch <- inter[posMismatcht]

    noMisIndividuals <- round(mismatchImplanted*implanted)


    for (i1 in posMismatch) {
      individind <- sample(implanted,noMisIndividuals)
      for (i1S in individind) {
        if ((i1>=intervalImp[i1S,3])&&(i1<=intervalImp[i1S,4])) {
          whichIndividual <- impl[i1S]
          if (alleleIimp[whichIndividual,i1]==0) {
            alleleIimp[whichIndividual,i1] <- 1
            alleleNimp[whichIndividual,i1] <- snvMinor[i1]
          } else {
            alleleIimp[whichIndividual,i1] <- 0
            alleleNimp[whichIndividual,i1] <- snvMajor[i1]
          }
        }

      }

    }
  }


}


alleleIimpGeno <- alleleIimp[s1,] + alleleIimp[s2,]

dummyL <- 1:snvs
dummy1 <- rep(1,snvs)
dummyC <- rep(0,snvs)
dummy <- as.character(dummyL)

annot <- list()

annot[[1]] <- dummy1
annot[[2]] <- pos
annot[[3]] <- snvNames
annot[[4]] <- snvMajor
annot[[5]] <- snvMinor
annot[[6]] <- dummy
annot[[7]] <- dummy
annot[[8]] <- dummy
annot[[9]] <- dummy
annot[[10]] <- freq
annot[[11]] <- as.numeric(dummyC)


names(annot) <- c("chromosome","position","snvNames","snvMajor","snvMinor","quality","pass","info","fields","frequency","changed")

annot <- as.data.frame(annot)

 save(snvs,haploN,individualsN,dataA,pos,pos1,snvNames,snvMajor,snvMinor,freq,annot,file=paste(fileprefix,run,"Anno.Rda",sep=""))

 save(impl,implG,inter,implantI,implantN,intervalImp,file=paste(fileprefix,run,"Impl.Rda",sep=""))

  # BEAGLE
  ########

  v1 <- rep(1,haploN)
  dim(v1) <- c(1,haploN)

  v0 <- as.numeric(rbind(dataA,dataA))
  dim(v0) <- c(1,haploN)



  dataB1 <- rbind(v0,v1,t(alleleNimp))

  col1 <- c("id","BC58",snvNames)

  dim(col1) <- c((snvs+2),1)

  col2 <- rep("M",snvs)

  col2 <- c("I","A",col2)

  dim(col2) <- c((snvs+2),1)


  dataB <- cbind(col2,col1,dataB1)

  write.table(dataB,file=paste(fileprefix,run,"beagle.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)




  # PLINK
  #######

  dataP1 <- matrix("character",nrow=individualsN,ncol=2*snvs)


  dataP1[,g1] <- alleleNimp[s1,]
  dataP1[,g2] <- alleleNimp[s2,]


  dataP <- cbind(plinkCols,dataP1)

  write.table(dataP,file=paste(fileprefix,run,"plink.ped",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)



  write.table(plinkmap,file=paste(fileprefix,run,"plink.map",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


  write.table(plinkCols,file=paste(fileprefix,run,"plink.fam",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")






   #MCMC
   #####

  write.table(alleleIimpGeno,file=paste(fileprefix,run,"mcmc.genotype",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  write.table(mcmc,file=paste(fileprefix,run,"mcmc.posmaf",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  write.table(initM,file=paste(fileprefix,run,"mcmc.initz",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


   #RELATE
   ########


  write.table((alleleIimpGeno+1),file=paste(fileprefix,run,"relate.geno",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

  write.table(pos1,file=paste(fileprefix,run,"relate.pos",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(chr,file=paste(fileprefix,run,"relate.chr",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)


  #FABIA
  ######


  namesL <-  cbind(1:haploN,1:haploN)
  write.table(namesL,file=paste(fileprefix,run,"fabia_individuals.txt",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)


  write(as.integer(haploN),file=paste(fileprefix,run,"fabia_annot.txt",sep=""),ncolumns=100)
  write(as.integer(snvs),file=paste(fileprefix,run,"fabia_annot.txt",sep=""),append=TRUE,ncolumns=100)
  write.table(annot,paste(fileprefix,run,"fabia_annot.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)


  write(as.integer(haploN),file=paste(fileprefix,run,"fabia_mat.txt",sep=""),ncolumns=100)
  write(as.integer(snvs),file=paste(fileprefix,run,"fabia_mat.txt",sep=""),append=TRUE,ncolumns=100)

  for (i in 1:haploN) {

    a1 <- which(alleleIimp[i,]>0.01)

    al <- length(a1)
    b1 <- alleleIimp[i,a1]

    a1 <- a1 - 1
    dim(a1) <- c(1,al)
    b1 <- format(as.double(b1),nsmall=1)
    dim(b1) <- c(1,al)

    write.table(al,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)
    write.table(a1,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)
    write.table(b1,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)

  }


}

}

split_sparse_matrix <- function(fileName,sparseMatrixPostfix="_mat.txt",intervalSize=10000,shiftSize=5000,annotation=TRUE) {

    narg <- as.integer(6)
    arg1 <- as.character(fileName)
    arg2 <- as.character(sparseMatrixPostfix)
    arg3 <- format(intervalSize, scientific = FALSE)
    arg4 <- format(shiftSize, scientific = FALSE)
    if (annotation) {
        arg5 <- as.character(1)
    } else {
        arg5 <- as.character(0)
    }

    message("Running 'split_sparse_matrix' on ",fileName)
    message("   Postfix of sparse matrix format ---- : ",sparseMatrixPostfix)
    message("   Interval size ---------------------- : ",intervalSize)
    message("   Shift size ------------------------- : ",shiftSize)
    if (annotation) {
    message("   Annotation is available ------------ : ")
    } else {
    message("   Annotation is not available -------- :")
    }

    .Call("split_sparse_matrix",narg,arg1,arg2,arg3,arg4,arg5,PACKAGE="hapFabia")

    message("")
    message("Split End.")


}


vcftoFABIA <- function(fileName,prefixPath="",noSnvs=NULL,outputFile=NULL) {


    arg1 <- as.character(fileName)
    arg2 <- as.character(prefixPath)
    if (is.null(noSnvs)) {
        arg3 <-  as.character("NA")
    } else {
        arg3 <- format(noSnvs, scientific = FALSE)
    }

    if (is.null(outputFile)) {
        arg4 <-  as.character("NA")
    } else {
        arg4 <- as.character(outputFile)
    }

    message("Running 'vcftoFABIA' on ",fileName)
    message("   Path to file ----------------------- : ",prefixPath)
    if (!is.null(noSnvs)) {
    message("   Number of SNVs --------------------- : ",noSnvs)
    } else {
    message("   Number of SNVs are unknown --------- :")
    }
    if (!is.null(outputFile)) {
    message("   Output file prefix------------------ : ",outputFile)
    } else {
    message("   Output file prefix given by input ----")
    }


    .Call("vcftoFABIA",arg1,arg2,arg3,arg4,PACKAGE="hapFabia")

    message("")
    message("Convert End.")


}

simulateIBDsegmentsFabia <- function(fileprefix="dataSim",minruns=1,maxruns=1,snvs=1000,individualsN=100,avDistSnvs=100,avDistMinor=10,noImplanted=1,implanted=10,length=50,minors=30,mismatches=0,mismatchImplanted=0.5,overlap=50) {


rateSnvs <- 1/avDistSnvs

rateMinor <- 1/avDistMinor



haploN <- 2*individualsN

dataA <- 1:individualsN



physPos <- vector("integer",snvs)


alleleI <- matrix(0,haploN,snvs)


#################

for (run in minruns:maxruns) {

#  write(noImplanted,file=paste(fileprefix,run,"Impl.txt",sep=""),ncolumns=100)


  distsSnvs <- rexp(snvs,rateSnvs)

  aa <- as.integer(0)
  for (i in 1:snvs) {
    aa <- aa + 1 + as.integer(round(distsSnvs[i]))
    physPos[i] <- aa

  }


  for (s in 1:haploN) {
    for (i in 1:snvs) {
      if (runif(1)<rateMinor) {
        alleleI[s,i] <- 1
      } else {
        alleleI[s,i] <- 0
      }
    }
  }



  freq <- colSums(alleleI)/haploN
  pos <- physPos
  pos1 <- pos/1000000


  alleleIimp <- alleleI

  allinter <- list()
  allimpl <- list()

  for (noIm in 1:noImplanted) {



  start <- sample((snvs-length-1),1)

  end <- start+length-1

  inter <- start:end

   posMinor <- sample(length,minors)

  implantI <- rep(0,length)

  implantI[posMinor] <- 1


  impl <- sample(haploN,implanted)

  allimpl[[noIm]] <- impl


  impl=sort(impl)

  implG <- (impl+1)%/%2
#  write(impl,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
#  write(implG,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
#  write(inter,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
#  write(implantI,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)


  intervalImp <- matrix(0,implanted,4)


  for (i in 1:implanted) {


    startI <- 1+sample((length-overlap)%/%2,1)
    endI <- length - 1 - sample((length-overlap)%/%2,1)
    interI <- inter[startI:endI]
    alleleIimp[impl[i],interI] <- implantI[startI:endI]

    interImp <- c(inter[startI],inter[endI],startI,endI)
    intervalImp[i,] <-  interImp

#    write(interImp,file=paste(fileprefix,run,"Impl.txt",sep=""),append=TRUE,ncolumns=100)
  }

  if (mismatches > 1) {

    posMismatcht <- 1+ sample((length-2),mismatches)

    posMismatch <- inter[posMismatcht]

    noMisIndividuals <- round(mismatchImplanted*implanted)


    for (i1 in posMismatch) {
      individind <- sample(implanted,noMisIndividuals)
      for (i1S in individind) {
        if ((i1>=intervalImp[i1S,3])&&(i1<=intervalImp[i1S,4])) {
          whichIndividual <- impl[i1S]
          if (alleleIimp[whichIndividual,i1]==0) {
            alleleIimp[whichIndividual,i1] <- 1
             } else {
            alleleIimp[whichIndividual,i1] <- 0
           }
        }

      }

    }
  }


}


dummy1 <- rep(1,snvs)
dummyL <- 1:snvs
dummyC <- rep(0,snvs)
dummy <- as.character(dummyL)
dummyM <- rep("A",snvs)
dummyU <- rep("T",snvs)
snvNames <- dummy

annot <- list()

annot[[1]] <- dummy1
annot[[2]] <- pos
annot[[3]] <- snvNames
annot[[4]] <- dummyM
annot[[5]] <- dummyU
annot[[6]] <- dummy
annot[[7]] <- dummy
annot[[8]] <- dummy
annot[[9]] <- dummy
annot[[10]] <- freq
annot[[11]] <- as.numeric(dummyC)


names(annot) <- c("chromosome","position","snvNames","snvMajor","snvMinor","quality","pass","info","fields","frequency","changed")

annot <- as.data.frame(annot)

namesL <-  cbind(1:haploN,1:haploN)

simu <- list(namesL,haploN,snvs,annot,alleleIimp)
save(simu,file=paste(fileprefix,run,"fabia_simulation.RData",sep=""))


write.table(namesL,file=paste(fileprefix,run,"fabia_individuals.txt",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)

write(as.integer(haploN),file=paste(fileprefix,run,"fabia_annot.txt",sep=""),ncolumns=100)
write(as.integer(snvs),file=paste(fileprefix,run,"fabia_annot.txt",sep=""),append=TRUE,ncolumns=100)
write.table(annot,paste(fileprefix,run,"fabia_annot.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)

# save(snvs,haploN,individualsN,dataA,pos,pos1,snvNames,snvMajor,snvMinor,freq,annot,file=paste(fileprefix,run,"Anno.Rda",sep=""))

# save(impl,implG,inter,implantI,implantN,intervalImp,file=paste(fileprefix,run,"Impl.Rda",sep=""))


  write(as.integer(haploN),file=paste(fileprefix,run,"fabia_mat.txt",sep=""),ncolumns=100)
  write(as.integer(snvs),file=paste(fileprefix,run,"fabia_mat.txt",sep=""),append=TRUE,ncolumns=100)

  for (i in 1:haploN) {

    a1 <- which(alleleIimp[i,]>0.01)

    al <- length(a1)
    b1 <- alleleIimp[i,a1]

    a1 <- a1 - 1
    dim(a1) <- c(1,al)
    b1 <- format(as.double(b1),nsmall=1)
    dim(b1) <- c(1,al)

    write.table(al,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)
    write.table(a1,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)
    write.table(b1,file=paste(fileprefix,run,"fabia_mat.txt",sep=""), sep = " ", quote = FALSE,row.names = FALSE,col.names = FALSE,append=TRUE)

  }


}

}

