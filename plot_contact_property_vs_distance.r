## Author: Qi Yu
## Contact: dryuqi@gmail.com


library(ggplot2)
library(HiTC)
library(scales)
library(dplyr)


gs<-function(genome){
  if (genome=="mm10"){
    genomelength<-read.table("/home/qiyu/data/genome/mm10/mouse.mm10_short.genome")}
  
  if (genome=="hg19"){
    genomelength<-read.table("/home/qiyu/data/genome/hg19/hg19.genome.new")
  }
  names(genomelength)<-c("chr","len")
  genomelength2<-data.frame(genomelength[,2],row.names=genomelength[,1])
  mylist<-list(genomelength,genomelength2)
  return(mylist)
}

plotIntraDistnew <- function(homo, winsize=NA, title=NA, genome="mm10",...){
  
  if (!is.na(winsize)){
    
    fnew<-data.frame(xp=double(),yp=double(),chr=character())
    k<-gs(genome)
    genomelength<-k[[1]]
    genomelength2<-k[[2]]
    
    for (ch in genomelength$chr){
      homog <- homo[which(homo$chr==ch),]
      if (nrow(homog)>0){
        r <- range(homog$distance)
        end<-round(logb(r[2]/winsize,1.12))
        if (end >0){
          #end<-122
          start<-r[1]
          test <- vector(mode="numeric", length=0)
          winn<-sapply(0:end,function(i){
            
            test[i+1]<-round(winsize*1.12^i)
            #return(test)
          })
          #fnew<-data.frame(xp=double(),yp=double(),chr=character())
          #win <- seq.int(from=r[1], to=r[2], by=winsize)
          xp <- yp <- rep(NA, length(winn)-1)
          
          tmp <- sapply(1L:(length(winn)-1L), function(i){
            idx <- which(homog$distance>=winn[i] & homog$distance<winn[i+1])
            #sidx<-which(homog$start>=winn[i] & homog$start<winn[i+1])
            #eidx<-which(homog$end>=winn[i] & homog$end<winn[i+1])
            #x <- mean(homog$distance[idx], na.rm=TRUE)
            x <- winn[i]
            #y<-0
            #if (length(homog$distance[idx])>5){
            # y <- length(homog$distance[idx])
            
            #/(genomelength2[c(ch),]-(winn[i+1]-winn[i])-1)
            #y <- length(homog$distance[idx])/(genomelength2[c(ch),]-winn[i]-1)
            y<-length(homog$distance[idx])/(length(homog$start)*(winn[i])/genomelength2[c(ch),])
            yb<-length(homog$distance[idx])
            #}
            return(c(x, y,yb))
          })
          
          xp <- tmp[1,]
          yp <- tmp[2,]
          ybp <- tmp[3,]
          
          names(xp) <- paste("n",1:length(xp), sep="")
          names(yp) <- paste("n",1:length(yp), sep="")
          names(ybp) <- paste("n",1:length(ybp), sep="")
          
          new<-as.data.frame(cbind(xp,yp,ybp))
          new<-new[which(new$yp>0),]
          if (length(new$xp)>0){
            new$chr<-ch
            fnew<-rbind(fnew,new)
          }#end length new >0
        }#end end>0
      }#end nrow >0 
    }#end per chr
  }#end if winsize is not NA
  write.table(fnew,row.names=F,file=paste(title,"_",genome,"_",winsize,"_fragment.out",sep=""))
  #write.table(yp,file="yp.out")
  doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -
                                                   min(x, na.rm=TRUE))}
  fnew$yp<-doit(fnew$yp)
  d=ggplot(fnew,aes(x=xp,y=yp,color=chr))+
    geom_line(aes(x=xp,y=yp,color=chr),size=1)+
    geom_point(size=1)+
    #geom_smooth()+
    xlab("Genomic Distance")+
    ylab("Contact Properties")+
    theme_bw()+
    theme(text = element_text(size=25))+
    ggtitle(title)+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides="lb")
  ggsave(plot=d,filename=paste(title,"_",genome,"_",winsize,"fragment_pos_per_chr.png", sep=''),height=10,width=16)
  
  fnew<-fnew[which(fnew$chr!="chrX" &fnew$chr!="chrY"),]
  fnew$yp<-doit(fnew$yp)
  d3=ggplot(fnew,aes(x=xp,y=yp,color=chr))+
    geom_line(size=1)+
    geom_point(aes(x=xp,y=yp,color=chr),size=1)+
    #geom_smooth()+
    xlab("Genomic Distance")+
    ylab("Contact Properties")+
    theme_bw()+
    theme(text = element_text(size=25))+
    ggtitle(title)+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides="lb")
  ggsave(plot=d3,filename=paste(title,"_",genome,"_",winsize,"fragment_pos_per_chr_noXY.png", sep=''),height=10,width=16)
  
  fnew1<-aggregate(fnew,by=list(fnew$xp),FUN=mean,na.rm=TRUE)
  write.table(fnew1,row.names=F,file=paste(title,"_",genome,"_",winsize,"_fragment_aggregate1.out",sep=""))
  fnew$yp<-doit(fnew$yp)
  d2=ggplot(fnew1,aes(x=xp,y=yp))+
    geom_line(aes(x=xp,y=yp),size=1)+
    geom_point(aes(x=xp,y=yp),size=1)+
    #geom_smooth()+
    xlab("Genomic Distance")+
    ylab("Contact Properties")+
    theme_bw()+
    theme(text = element_text(size=25))+
    ggtitle(title)+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
    annotation_logticks(sides="lb")
  ggsave(plot=d2,filename=paste(title,"_",genome,"_",winsize,"fragment_pos1.png", sep=''),height=10,width=10)
  
}#end of function plotIntraDistnew

readfile<-function(filename,title=NA,genome="mm10"){
  homo<-read.table(filename)
  names(homo)<-c("chr","start","end")
  homo$distance<-homo$end-homo$start
  plotIntraDistnew(homo,winsize=40000,title=title,genome=genome)
}#end of function readfile

#readfile("/home/yuqi/Downloads/GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5_hg19.validPair.txt.gz_count.txt",title="hela_test",genome="hg19")
#readfile("/home/yuqi/Downloads/GSM2745897_18FEB15_PE50_C66B1AC-A_Sample_HiC1-aka-CCHiC-HeLaS3CCL2p2-M-98a_hg19.validPair.txt.gz_count.txt",title="Hela_new",genome="hg19")
#readfile("/mnt/RDCO/shared/Gang/sperm/aligned/merged_nodups.txt_count.txt",title="sperm",genome="mm10")
#readfile("/mnt/RDCO/shared/Gang/aligned_Lep/merged_nodups.txt_count.txt",title="Lep",genome="mm10")
#readfile("/mnt/RDCO/shared/Gang/aligned_preLep/merged_nodups.txt_count.txt",title="preLep",genome="mm10")
#readfile("/mnt/RDCO/shared/Gang/aligned_4N_ZPD/merged_nodups.txt_count.txt",title="4N_ZPD",genome="mm10")

readfile("/mnt/RDCO/shared/Gang//aligned_spermatid/merged_nodups.txt_count.txt",title="spermatid",genome="mm10")
readfile("/mnt/RDCO/shared/Gang//aligned_zygotene/merged_nodups.txt_count.txt",title="zygotene",genome="mm10")
readfile("/mnt/RDCO/shared/Gang//aligned_pachytene/merged_nodups.txt_count.txt",title="pachytene",genome="mm10")
readfile("/mnt/RDCO/shared/Gang//aligned_diplotene/merged_nodups.txt_count.txt",title="diplotene",genome="mm10")
#readfile("./HiCPro/hicpro_test/hic_results/data/dixon_2M_2/SRR400264_01_hg19.bwt2pairs.validPairs.bed",title="test_human",genome="hg19")

#l <- sapply(list.files("/home/qiyu/data/hiC/lep", pattern=paste("inter_30_"), full.names=TRUE),import.my5C)
#hiC <- HTClist(l)
#CQC1(hiC, winsize = 1e+06, dev.new=FALSE, hist.dist=FALSE,hist.interac=FALSE)
