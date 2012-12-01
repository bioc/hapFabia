### ---------------------------
### Factorization class methods
### ---------------------------

setMethod("topLZ",signature(res="Factorization",n="numeric",LZ="character",indices="logical",p="ANY",w="ANY"),
function(res,n=1,LZ="L",indices=TRUE,p=NULL,w=NULL) {

if (missing(res)) {
    stop("Fabia result 'res' is missing. Stopped.")
}

if ((is.null(p))&&(is.null(w))) {

    stop("Either 'p' or 'w' must be given. Stopped.")

}

if (LZ=="Z") {

    if (!is.null(p)) {
        if (indices) {
            return(which(Z(res)[n,]>quantile(Z(res)[n,],p)))
        } else {
          return(Z(res)[n,which(Z(res)[n,]>quantile(Z(res)[n,],p))])
        }
    } else {
        if (indices) {
            return(which(Z(res)[n,]>w))
        } else {
            return(Z(res)[n,which(Z(res)[n,]>w)])
        }

    }

} else {

    if (!is.null(p)) {
        if (indices) {
            return(which(L(res)[,n]>quantile(L(res)[,n],p)))
        } else {
          return(L(res)[which(L(res)[,n]>quantile(L(res)[,n],p)),n])
        }
    } else {
        if (indices) {
            return(which(L(res)[,n]>w))
        } else {
            return(L(res)[which(L(res)[,n]>w),n])
        }

    }
}


}
)


setMethod("plotL",signature(res="Factorization",n="numeric",p="ANY",w="ANY",type="character",intervv="numeric",off="numeric",t="character",cex="numeric"),
function(res,n=1,p=NULL,w=NULL,type="points",intervv=500,off=0,t="p",cex=1) {

if (missing(res)) {
    stop("Fabia result 'res' is missing. Stopped.")
}

if ((is.null(p))&&(is.null(w))) {

    stop("Either 'p' or 'w' must be given. Stopped.")

}

if (!is.null(p)) {

    if (type=="points") {
        a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=p,w=NULL)
        b <- topLZ(res=res,n=n,LZ="L",indices=FALSE,p=p,w=NULL)
        plot(a,b,xlim=c(0,dim(L(res))[1]),col="blue",type=t,cex=cex,main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

    } else {
        if (type=="histogram") {
            off <- off%%intervv
            bb <- seq(1-off,length(L(res)[,n]),intervv)
            lb <- length(bb)
            bb <- c(bb,(bb[lb]+intervv))
            a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=p,w=NULL)
            aa <- hist(a,breaks=bb,plot = FALSE)
            plot(aa,freq = TRUE,xlim=c(0,dim(L(res))[1]),col="blue",main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

        } else {
            a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=p,w=NULL)
            b <- topLZ(res=res,n=n,LZ="L",indices=FALSE,p=p,w=NULL)
            smoothScatter(a,b,xlim=c(0,dim(L(res))[1]),main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

        }

    }


} else {

    if (type=="points") {
        a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=NULL,w=w)
        b <- topLZ(res=res,n=n,LZ="L",indices=FALSE,p=NULL,w=w)
        plot(a,b,xlim=c(0,dim(L(res))[1]),col="blue",type=t,cex=cex,main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

    } else {
        if (type=="histogram") {
            off <- off%%intervv
            bb <- seq(1-off,length(L(res)[,n]),intervv)
            lb <- length(bb)
            bb <- c(bb,(bb[lb]+intervv))
            a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=NULL,w=w)
            aa <- hist(a,breaks=bb,plot = FALSE)
            plot(aa,freq = TRUE,xlim=c(0,dim(L(res))[1]),col="blue",main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

        } else {
            a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=NULL,w=w)
            b <- topLZ(res=res,n=n,LZ="L",indices=FALSE,p=NULL,w=w)
            smoothScatter(a,b,xlim=c(0,dim(L(res))[1]),main=paste("bicluster ",n,sep=""),ylab="counts",xlab="SNV number")

        }

    }


}


}
)

setMethod("histL",signature(res="Factorization",n="numeric",p="ANY",w="ANY",intervv="numeric",off="numeric"),
function(res,n=1,p=NULL,w=NULL,intervv=500,off=0) {

    if (missing(res)) {
        stop("Fabia result 'res' is missing. Stopped.")
    }

    if ((is.null(p))&&(is.null(w))) {

        stop("Either 'p' or 'w' must be given. Stopped.")

    }

    off <- off%%intervv
    bb <- seq(1-off,length(L(res)[,n]),intervv)
    lb <- length(bb)
    bb <- c(bb,(bb[lb]+intervv))
    if (!is.null(p)) {
        a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,p=p,w=NULL)
    } else {
        a <- topLZ(res=res,n=n,LZ="L",indices=TRUE,w=w,p=NULL)
    }
    rhist <- hist(a,breaks=bb,plot = FALSE)
    return(rhist)
}
)





setMethod("extractIBDsegments",signature(res="Factorization",sPF="list",annot="data.frame",chrom="character",labelsA="matrix",ps="numeric",psZ="numeric",inteA="numeric",thresA="numeric",mintagSNVs="numeric",off="numeric",procMinIndivids="numeric",thresPrune="numeric"),
function(res,sPF,annot=NULL,chrom="",labelsA=NULL,ps=0.9,psZ=0.8,inteA=500,thresA=11,mintagSNVs=8,off=0,procMinIndivids=0.1,thresPrune=1e-3) {

    if (missing(res)) {
        stop("Fabia result 'res' is missing. Stopped.")
    }

    if (missing(sPF)) {
        stop("Individuals per SNV (genotype data) 'sPF' is missing. Stopped.")
    }

if (is.null(labelsA)) {
  nnL <- length(Z(res)[1,])
  labelsA <- cbind(as.character(1:nnL),as.character(1:nnL),as.character(1:nnL),as.character(1:nnL))
}


if (is.null(annot)) {

  snvs <- length(L(res)[,1])

  dummyL <- 1:snvs
  dummyC <- rep(0,snvs)
  dummyN <- rep("N",snvs)
  dummy <- as.character(dummyL)

  annot <- list()

  annot[[1]] <- dummyL
  annot[[2]] <- 100*dummyL
  annot[[3]] <- dummy
  annot[[4]] <- dummyN
  annot[[5]] <- dummyN
  annot[[6]] <- dummy
  annot[[7]] <- dummy
  annot[[8]] <- dummy
  annot[[9]] <- dummy
  annot[[10]] <- as.double(dummyC)
  annot[[11]] <- dummyC

}

IBDsegment_res <- new("IBDsegmentList",IBDsegments=list(),lengthList=0)
IBDsegment_res_idx <- 0

max_n <- ncol(L(res))

for (n in 1:max_n) {

if (max(L(res)[,n])>0.000001) {


ib <- findDenseRegions(L(res)[,n],p=ps,inte=inteA,thres=thresA,off=off)

topZ <- which(Z(res)[n,]>quantile(Z(res)[n,],psZ))


if ((!is.null(ib$len))&&(length(ib$len)>0)) {

for (i in 1:length(ib$len)) {


ibb <- as.vector(unlist(ib$pos[[i]]))
#ma <- length(ibb)

il <- length(ibb)


aa <- c()
nu <- list()
for (j in 1:il) {
  nu[[j]] <- sPF$sL[[ibb[j]]]
  aa <- c(aa,nu[[j]])
}

aa <- unique(aa)

aa <- intersect(aa,topZ)

lq <- length(aa)

if (lq>1) {

aa <- sort(aa)



mat <- matrix(0,lq,il)

for (j in 1:il) {
  nv <- which(match(aa,nu[[j]])>0)
  mat[nv,j] <- rep(1,length(nv))
}



select1 <- which(rowSums(mat)>mintagSNVs)

lq1 <- length(select1)


if (lq1>1) {

individual1 <- aa[select1]



mat1 <- mat[select1,]

select2 <- which(colSums(mat1)>1)

lq2 <- length(select2)

if (lq2>1) {

mat1 <- mat1[,select2]


tagSNV1 <- ibb[select2]



matB <- mat1
ig=1
ie=1
while (ie == 1) {


    nnc <- ncol(matB)
    subclS <- 1:nnc
    lq2A <- nnc


    nnl <- nrow(matB)
    subclM <- 1:nnl
    lq1A <- nnl



    lq2A_old <- lq2A
    lq1A_old <- lq1A


    done=0

    if ((lq1A<2)||(lq2A<2)) {
      lq1A <- 1
      done <- 1
    }

    while (done==0) {




      matA <- matB[subclM,subclS]

      cs1 <- rowSums(matA)

      m2 <- which.max(cs1)

      v1 <- matA%*%matA[m2,]

      v1A <- v1
      v1A[m2] <- 0

      m3 <- which.max(v1A)

      v2 <-  matA%*%matA[m3,]

      subclM1A <- which(v1>mintagSNVs)
      subclM1B <- which(v2>mintagSNVs)

      subclM1 <- intersect(subclM1A,subclM1B)

      lq1A <- length(subclM1)



     if (lq1A>1)  {



      cs2 <- colSums(matA[subclM1,])

      subclS1 <- which(cs2>max(1,(lq1A%/%(1/procMinIndivids))))

      lq2A <- length(subclS1)



      subclM <- subclM[subclM1]
      subclS <- subclS[subclS1]

      subclM2 <- subclM[c(m2,m3)]




      if ((lq1A_old==lq1A)&&(lq2A_old==lq2A)) {
        done <- 1
      } else {
        lq2A_old <-  lq2A
        lq1A_old <-  lq1A
      }

       if ((lq1A<2)||(lq2A<2)) {
        lq1A <- 1
        done <- 1
      }

    } else {
       done <- 1
    }

    }


    if (lq1A>1) {


      a0 <- tagSNV1[subclS]
      a1 <- diff(a0)
      m <- 1/(max(2,median(a1)))
      k0 <- length(a1)
      k <- k0
      wegS <- 0
      wegE <- 0
      thresPruneN <- 1.0 - thresPrune
      if (k>2) {
        doPrune <- TRUE
      } else {
        doPrune <- FALSE
      }

      while (doPrune) {

        doPrune <- FALSE

        if ((k>2)&&(pexp(a1[1],m) > thresPruneN)) {
          a1 <- a1[-1]
          k <- k - 1
          m <- 1/(max(2,median(a1)))
          wegS <- wegS+1
          doPrune <- TRUE
        } else {
          if ((k>2)&&(pexp(a1[2],m) > thresPruneN)) {
            a1 <- a1[c(-2,-1)]
            k <- k - 2
            m <- 1/(max(2,median(a1)))
            wegS <- wegS+2
            doPrune <- TRUE
          }

        }

        if ((k>2)&&(pexp(a1[k],m) > thresPruneN)) {
          a1 <- a1[-k]
          k <- k - 1
          m <- 1/(max(2,median(a1)))
          wegE <- wegE+1
          doPrune <- TRUE
        } else {
          if ((k>2)&&(pexp(a1[k-1],m) > thresPruneN)) {
            a1 <- a1[c(-k,(-k+1))]
            k <- k - 2
            m <- 1/(max(2,median(a1)))
            wegE <- wegE+2
            doPrune <- TRUE
          }

        }

      }

      subclS <- subclS[(1+wegS):(k0-wegE)]
      lq2A <- length(subclS)



      if (lq2A>=mintagSNVs) {


      matA <- matB[subclM,subclS]


      individual <- as.vector(individual1[subclM])
      individual2 <- as.vector(individual1[subclM2])

      labelsEUA <- as.vector(labelsA[individual,2])
      labelsNAA <- as.vector(labelsA[individual,1])
      labelsPlatformA <- as.vector(labelsA[individual,4])
      labels_ALLA <- as.vector(paste(labelsNAA,labelsEUA,sep="_"))

      labelsNA2 <- as.vector(labelsA[individual2,1])

      tagSNV <- as.vector(tagSNV1[subclS])


      physRangeA <- as.vector(annot[[2]][tagSNV])
      IBDsegmentLength <-  physRangeA[lq2A]-physRangeA[1]
      IBDsegmentPos <- (physRangeA[lq2A]+physRangeA[1])%/%2
      allelesA <- as.vector(paste(annot[[4]][tagSNV],":",annot[[5]][tagSNV],sep=""))
      tagSNVNamesA <- as.vector(annot[[3]][tagSNV])
      tagSNVFreqA <- as.vector(annot[[10]][tagSNV])
      nnn <- length(labelsA[,1])
      tagSNVGroupFreqA <- as.vector(sPF$nsL[tagSNV]/nnn)
      tagSNVChangeA <- as.vector(annot[[11]][tagSNV])

      tagSNVsPerIndividual <- as.vector(rowSums(matA))
      individualPerTagSNV <- as.vector(colSums(matA))
      tagSNVAnnoA <- rep(as.list(""),length(tagSNV))

      IBDsegment_res_idx = IBDsegment_res_idx + 1


      chrom <- names(which.max(table(annot[[1]][tagSNV])))[1]

      if ( !is.character(chrom) || (nchar(chrom)==0) || (nchar(chrom)>2 )) {
          chrom <- ""
      }

      res_t <- new("IBDsegment",ID=IBDsegment_res_idx,bicluster_id=n,chromosome=chrom,IBDsegmentPos=IBDsegmentPos,IBDsegmentLength=IBDsegmentLength,numberIndividuals=lq1A,numbertagSNVs=lq2A,individuals=individual,tagSNVs=tagSNV,populationIndividuals=labelsEUA,idIndividuals=labelsNAA,labelIndividuals=labels_ALLA,platformIndividuals=labelsPlatformA,coreClusterIndividuals=labelsNA2,tagSNVPositions=physRangeA,tagSNVAlleles=allelesA,tagSNVNames=tagSNVNamesA,tagSNVFreq=tagSNVFreqA,tagSNVGroupFreq=tagSNVGroupFreqA,tagSNVChange=tagSNVChangeA,tagSNVsPerIndividual=tagSNVsPerIndividual,individualPerTagSNV=individualPerTagSNV,tagSNVAnno=tagSNVAnnoA)

      IBDsegment_res[[IBDsegment_res_idx]] <- res_t


        ig <- ig + 1
        if ((lq1A<(nnl-1))&&(lq2A<(nnc-1))) {

            matB <-  matB[-subclM,-subclS]

            individual1 <- individual1[-subclM]


            tagSNV1 <- tagSNV1[-subclS]

       } else {
            ie=100
        }

        } else {
            ie=100
        }


    } else {
        ie=100
    }




  }

}

}

}
}
}
}

}


return(IBDsegment_res)

}
)

