## Author: Qi Yu
## Contact: dryuqi@gmail.com


library(ggplot2)
library(HiTC)

plotIntraDist <- function(xdata.intra, title=title,xdata.intra.dist, trim.range=.98, winsize=NA, add=FALSE, log=TRUE, fit=FALSE, fit.lim=NA, fit.out=1, ...){
  stopifnot(length(xdata.intra)==length(xdata.intra.dist))
  r <- range(xdata.intra.dist)
  
  ## Remove zeros
  idx <- which(xdata.intra>0)
  xdata.intra <- xdata.intra[idx]
  xdata.intra.dist <- xdata.intra.dist[idx]
  
  ## Windowing
  ## To optimize - see tapply function
  if (!is.na(winsize)){
    
    win <- seq.int(from=r[1], to=r[2], by=winsize)
    xp <- yp <- rep(NA, length(win)-1)
    
    tmp <- sapply(1L:(length(win)-1L), function(i){
      idx <- which(xdata.intra.dist>=win[i] & xdata.intra.dist<win[i+1])
      
      sub <- xdata.intra.dist[idx]
      x <- mean(xdata.intra.dist[idx], na.rm=TRUE)
      y <- mean(xdata.intra[idx], na.rm=TRUE)
      return(c(x, y))
    })
    
    xp <- tmp[1,]
    yp <- tmp[2,]
    ptype <- "b"
  }else{
    xp <- xdata.intra.dist
    yp <- xdata.intra
    ptype <- "p"
  }
  names(xp) <- paste("n",1:length(xp), sep="")
  names(yp) <- paste("n",1:length(yp), sep="")
  
  ## Trim the interaction counts
  if (trim.range<1){
    qt <- quantile(yp, probs=c((1-trim.range),trim.range), na.rm=TRUE)
    idx <- which(yp>=qt[1] & yp<=qt[2])
    yp <- yp[idx]
    xp <- xp[idx]
  }
  
  ## Log
  #if (log){
   # xp <- log10(xp)
  #  yp <- log10(yp)
  #}
  
  if (fit){
    ## remove extreme distances
    #qt <- quantile(xp, probs=c(0.1,0.8), na.rm=TRUE)
    #idx <- which(xp>=qt[1] & xp<=qt[2])
    #xp <- xp[idx]
    #yp <- yp[idx]
    
    ## Define the scaling region
    if (length(fit.lim)==2){
      yp.fit <- yp[which(xp>fit.lim[1] & xp<fit.lim[2])]
      xp.fit <- xp[which(xp>fit.lim[1] & xp<fit.lim[2])]
    }else{
      xp.fit <- xp
      yp.fit <- yp
    }
    res.fit <- lm(yp.fit~xp.fit)
    
    ## remove outliers on the fitted region
    ## outliers are defined as the points with the higher residual values (distance to the regression curve)
    if (fit.out<1){
      r <- residuals(res.fit)
      th<-quantile(r, probs=fit.out)
      out <- names(which(r>th))
      
      xp.fit <- xp.fit[setdiff(names(xp.fit), out)]
      yp.fit <- yp.fit[setdiff(names(xp.fit), out)]
      res.fit <- lm(yp.fit~xp.fit)
      
      xp <- xp[setdiff(names(xp), out)]
      yp <- yp[setdiff(names(xp), out)]
    }
  }
 
  ## Plotting function
  if (!add){
    new<-as.data.frame(cbind(xp,yp))
    #rownames
    write.table(new,row.names=F,file=paste(title,"c.out",sep=""))
    #write.table(yp,file="yp.out")
    d=ggplot(new)+
      geom_line(aes(x=new$xp,y=new$yp),size=2)+
      geom_point(aes(x=new$xp,y=new$yp),size=2)+
      xlab("Genomic Distance")+
      ylab("Interaction Counts")+
      theme_bw()+
      theme(text = element_text(size=25))+
      ggtitle(title)+
      scale_x_log10()+scale_y_log10()+
      annotation_logticks(sides="lb")
    ggsave(plot=d,filename=paste(title,"_1M.png", sep=''),height=6,width=6)

    plot(x=c(min(xp), max(xp)), y=c(min(yp), max(yp)),  xlab="Genomic Distance (log10)", ylab="Interaction Counts (log10)",
         frame=FALSE, type="n", ...)
  }
  
  if (fit){
    pcol <- RColorBrewer::brewer.pal(8, "Pastel2")
    if (length(fit.lim)==2)
      rect(fit.lim[1], min(yp)-1, fit.lim[2], max(yp), col=pcol[5], border=pcol[5])
    abline(res.fit, ...)
    #text(x=max(xp), y=max(yp), labels=paste("a=", round(res.fit$coefficients[2],6)), font=2, cex=.7)
  }
  points(x=xp, y=yp, type=ptype, ...)
  if (fit)
    return(res.fit)
  else
    invisible(NULL)
}



extractCounts <- function(x){
  x.intra <- x.inter <- NULL
  xdata.intra <- xdata.inter <- NULL
  
  if (inherits(x,"HTClist")){
    ## Separate Intra/Inter chromosomal interaction
    if (length(which(isIntraChrom(x)))>0){
      x.intra <- x[isIntraChrom(x)]
      xdata.intra <- unlist(lapply(x.intra,function(x){
        return(as(as(intdata(x), "sparseMatrix"),"dgTMatrix")@x)
      }))
      xdata.intra.dist <- unlist(lapply(x.intra,function(x){
        return(intervalsDist(x, use.zero=FALSE)@x)}))
    }
    if (length(which(!isIntraChrom(x)))>0){
      x.inter <- x[!isIntraChrom(x)]
      xdata.inter <- unlist(lapply(x.inter,function(x){
        return(intdata(x)@x)
      }))
    }
  } else if(inherits(x,"HTCexp")){
    if (isIntraChrom(x)){
      x.intra <- x
      xdata.intra <- as(as(intdata(x.intra), "sparseMatrix"),"dgTMatrix")@x
      xdata.intra.dist <- intervalsDist(x.intra, use.zero=FALSE)@x
    }else {
      x.inter <- x
      xdata.inter <- intdata(x.inter)@x
      xdata.intra.dist <- NULL
    }
  }else{
    stop("Wrong input type. 'HTCexp' or 'HTClist' objects expected.")
  }
  return(list(intra=xdata.intra, intra.dist=xdata.intra.dist, inter=xdata.inter))
}


###################################
## CQC
##
## 'C' Quality Control
##
## x = an object of class HTCexp/HTClist. In case of list, the x are merged as one.
## trans.ratio = if true, plot the histogram of inter/intrachromosomal interactions 
## hist.interac = if true, plot the interaction frequency. How many interaction have n reads 
## scat.interac.dist = if true, plot the scatter plot of counts vs genomic distances. Interaction Distance vs interaction frequency
## hist.dist = if true, plot an histogram of distances between y and x intervals. How many interaction have n distance
## trim.range = remove the extreme values by trimming the counts. Only used for plotting functions
## dev.new = if true, draw each plot in a new view
##
###################################

CQC1<- function(x, title,cis.trans.ratio=FALSE, hist.interac=FALSE, scat.interac.dist=TRUE, hist.dist=FALSE, trim.range=0.98, winsize=NA, dev.new=TRUE){
  
  message("Get data ...")
  data <- extractCounts(x)
  xdata.inter <- data$inter
  xdata.intra <- data$intra
  xdata <- c(xdata.inter, xdata.intra)
  xdata.intra.dist <- data$intra.dist
  
  message("Generate quality control plots ...")
  nbplot <- length(which(c(cis.trans.ratio, scat.interac.dist, hist.dist, hist.interac)))
  if (!is.null(xdata.inter)){
    nbplot <- nbplot+1
  }
  if (!dev.new)
    par(mfrow=c(ceiling(nbplot/2),2), mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
  

  ## Scatter plot of interaction distance for intrachromosomal interactions
  if(scat.interac.dist && !is.null(xdata.intra)>0){
    if (dev.new){
      dev.new()
      par(mar=c(4.1, 4.1, 2.5, 1.5), font.lab=2)
    }
    plotIntraDist(xdata.intra, title=title,xdata.intra.dist, trim.range, winsize=winsize,
                  cex=0.5, cex.lab=0.7, pch=20,
                  cex.axis=0.7, cex.main=0.9, main="Scatter Plot (Frequency(Y) vs Distance(X))\nCIS Interaction Counts")
  }
  
}

#exDir <- system.file("extdata", package="HiTC")
homo<-importC(con="/home/qiyu/data/hiC/HiCPro/allele_spe/hic_results/data/171031_HiC_4N_B6xCAST/171031_HiC_4N_B6xCAST_500000.matrix", 
        xgi="/home/qiyu/data/hiC/HiCPro/allele_spe/hic_results/data/171031_HiC_4N_B6xCAST/171031_HiC_4N_B6xCAST_500000_abs.bed")
khic<-HTClist(homo)

CQC1(khic, title="Homologous_pair",winsize = 1e+06, dev.new=FALSE, hist.dist=FALSE,hist.interac=FALSE)

#l <- sapply(list.files("/home/qiyu/data/hiC/lep", pattern=paste("inter_30_"), full.names=TRUE),import.my5C)
#hiC <- HTClist(l)
#CQC1(hiC, winsize = 1e+06, dev.new=FALSE, hist.dist=FALSE,hist.interac=FALSE)
