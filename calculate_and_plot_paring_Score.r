## Author: Qi Yu
## Contact: dryuqi@gmail.com


library("KernSmooth")
library(raster)
library(ggplot2)
library(tibble)
library(reshape2)
library(plyr)


convertm2<-function(file){
  
  plots1<-read.table(file,sep=" ",header=FALSE)
  names(plots1)<-c("chra","starta","enda","stranda","chrb","startb","endb","strandb")
  heto<-plots1[which(plots1$chra!=plots1$chrb),]
  homo<-plots1[which(plots1$chra==plots1$chrb),]
  newmean<-function(i){
  pos=with(i,cbind(starta,startb))
  chrlen = max(pos)
  gridsize = ceiling(chrlen/5e5)
  bandwidth = 5e5
  den = bkde2D(pos, bandwidth=c(1,1)*bandwidth, gridsize=c(1,1)*gridsize)
  den$fhat <- den$fhat + t(den$fhat)
  r <- as(den$fhat, "RasterLayer")
  rmean <- focal(r, w=matrix(1/49, ncol=7, nrow=7), fun=mean)
  meand<-diag(as(rmean,"matrix"))}
  
  newmeanhe<-newmean(heto)
  newmeanho<-newmean(homo)
  return(list(newmeanhe,newmeanho))
}#end of function readfile

homomean<-convertm2("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1.txt")
#50000

homomean2<-convertm2("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1.txt")
#2000000
homomean3<-convertm2("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1.txt")

convertm_divid<-function(file){
  
  plots1<-read.table(file,sep=" ",header=FALSE)
  names(plots1)<-c("chra","starta","enda","stranda","chrb","startb","endb","strandb")
  heto<-plots1[which(plots1$chra!=plots1$chrb),]
  homo<-plots1[which(plots1$chra==plots1$chrb),]
  
  newmean<-function(i){
    pos=with(i,cbind(starta,startb))
    chrlen = max(pos)
    gridsize = ceiling(chrlen/5e5)
    bandwidth = 5e5
    den = bkde2D(pos, bandwidth=c(1,1)*bandwidth, gridsize=c(1,1)*gridsize)
    den$fhat <- den$fhat + t(den$fhat)
    r <- as(den$fhat, "RasterLayer")
    rmean <- focal(r, w=matrix(1/49, ncol=7, nrow=7), fun=mean)
    meand<-diag(as(r,"matrix"))/diag(as(rmean,"matrix"))}
  
  newmeanhe<-newmean(heto)
  newmeanho<-newmean(homo)
  return(list(newmeanhe,newmeanho))
}#end of function readfile

homomean_d<-convertm_divid("181016_9014_3_CTTGTA_GC_050218_Hi_C_B6XCAST_Preleptotene.mm10_genome_mask_CAST.bwt2pairs_VI_chr1.txt")
homomean<-homomean_d

heto<-data.frame(log(homomean[[1]]))
names(heto)<-c("value")
heto$type<-"heto"
heto<-rownames_to_column(heto,var="location")
heto$location<-strtoi(heto$location)*500

homo<-data.frame(log(homomean[[2]]))
names(homo)<-c("value")
homo$type<-"homo"
homo<-rownames_to_column(homo,var="location")
homo$location<-strtoi(homo$location)*500

#heto$value<-log2(heto$value-min(heto$value))
#homo$value<-log2(homo$value-min(homo$value))


new<-rbind(heto,homo)
new<-na.omit(new)


newd<-dcast(new,location~type)
names(newd)<-c("end",'heto','homo')

EI<-read.table("prelep_B6XCAST.hic_EI.bedgraph_chr1")
names(EI)<-c("chr","start","end","EI")
EIs<-EI[which(EI$chr=="chr1"),]
EIs$start<-EIs$start/1000
EIs$end<-EIs$end/1000

test=data.frame(join_all(list(newd,EIs),by=c("end")))
test<-na.omit(test)
newtest<-melt(test[c(1,2,3,6)],id.vars="end")

d2=ggplot(newtest)+
  #geom_density(linesize=2)+
  geom_point(aes(x=end,y=value),size=2)+
  geom_line(aes(x=end,y=value),size=1)+
  theme_bw()+
  facet_wrap( ~ variable,ncol=1,scales = "free")+
  xlab("distance (kb)")+
  #xlim(20000,60000)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("chr1_heto_homo_EI_end_divid.png", sep=''),height=16,width=26)

d2=ggplot(newtest)+
  #geom_density(linesize=2)+
  geom_point(aes(x=end,y=value,color=variable),size=2)+
  geom_line(aes(x=end,y=value,color=variable),size=1)+
  theme_bw()+
  #facet_wrap( ~ variable,ncol=1,scales = "free")+
  xlab("distance (kb)")+
  xlim(40000,80000)+
  #ylim(-1,3)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("chr1_heto_homo_EI_end_divid_differcolor.png", sep=''),height=10,width=10)

d2=ggplot(test)+
  #geom_density(linesize=2)+
  geom_point(aes(x=heto,y=EI),size=2)+
  #geom_line(aes(x=start,y=value),size=1)+
  theme_bw()+
  #facet_wrap( ~ variable,ncol=1,scales = "free")+
  #xlab("distance (kb)")+
  #xlim(20000,60000)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

ggsave(plot=d2,filename=paste("chr1_heto_EI_end_divid.png", sep=''),height=6,width=8)

d2=ggplot(new)+
  #geom_density(linesize=2)+
  geom_point(aes(x=location,y=value,color=type),size=2)+
  geom_line(aes(x=location,y=value,color=type,group=type),size=1)+
  theme_bw()+
  #facet_wrap( ~ scale,ncol=5)+
  xlab("distance (kb)")+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

# axis.title.x=element_blank(),
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank(),
# axis.title.y=element_blank(),
# axis.text.y=element_blank(),
# axis.ticks.y=element_blank(),
#legend.title=element_blank(),legend.justification=c(1,01)

#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("Pairingscore_prelep_highres_scale.png", sep=''),height=16,width=26)

homomean_prel<-convertm2("180530_prelep_B6XCAST_genome_mask_CAST.bwt2pairs_homozygous.bam_ncc.txt_chr1.txt")


heto<-data.frame(homomean_prel[[1]])
names(heto)<-c("value")
heto$type<-"heto"
heto<-rownames_to_column(heto,var="location")
heto$location<-strtoi(heto$location)*400

homo<-data.frame(homomean_prel[[2]])
names(homo)<-c("value")
homo$type<-"homo"
homo<-rownames_to_column(homo,var="location")
homo$location<-strtoi(homo$location)*400

new<-rbind(heto,homo)
new<-na.omit(new)



rt1<-read.table("log2.ratio.preL_2N.mm10_rep_timing.bedgraph")
names(rt1)<-c("chr","start","end","timing")
rt1s<-rt1[which(rt1$chr=="chr1"),]
rt1s$start<-rt1s$start/1000

d2=ggplot(rt1s)+
  #geom_density(linesize=2)+
  geom_point(aes(x=start,y=timing),size=2)+
  geom_line(aes(x=start,y=timing),size=1)+
  theme_bw()+
  #facet_wrap( ~ scale,ncol=5)+
  xlab("distance (kb)")+
  xlim(20000,60000)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("replication_timing_scale3.png", sep=''),height=16,width=26)

d2=ggplot(new)+
  #geom_density(linesize=2)+
  geom_point(aes(x=location,y=value,color=type),size=2)+
  geom_line(aes(x=location,y=value,color=type,group=type),size=1)+
  theme_bw()+
  #facet_wrap( ~ scale,ncol=5)+
  xlab("distance (kb)")+
  xlim(20000,60000)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

# axis.title.x=element_blank(),
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank(),
# axis.title.y=element_blank(),
# axis.text.y=element_blank(),
# axis.ticks.y=element_blank(),
#legend.title=element_blank(),legend.justification=c(1,01)


#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("Pairingscore_prelep_highres_scale3.png", sep=''),height=16,width=26)

new<-merge(heto,homo,by="location")
new$ratio<-log(new$value.x/new$value.y)
new<-na.omit(new)


d2=ggplot(new)+
  #geom_density(linesize=2)+
  geom_point(aes(x=location,y=ratio),size=2)+
  geom_line(aes(x=location,y=ratio),size=1)+
  #geom_line(aes(x=rt1s$start,y=rt1s$timing),size=1)+
  theme_bw()+
  #facet_wrap( ~ scale,ncol=5)+
  xlab("distance (kb)")+
  xlim(20000,60000)+
  #scale_y_log10()+
  #annotation_logticks(sides="b")+
  theme(text = element_text(size=25))

# axis.title.x=element_blank(),
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank(),
# axis.title.y=element_blank(),
# axis.text.y=element_blank(),
# axis.ticks.y=element_blank(),
#legend.title=element_blank(),legend.justification=c(1,01)


#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename=paste("Pairingscore_prelep_highres_scale_ratio.png", sep=''),height=16,width=26)



plots1<-read.table("180901_9011_1_ACTGAT_GC_072318_Hi-C_B6XCAST_Zygotene.mm10_genome_mask_CAST.bwt2pairs_heterozygous.bam_ncc.txt_chr1.out.count.txt",sep=" ",header=FALSE)
names(plots1)<-c("chra","starta","enda","stranda","chrb","startb","endb","strandb")
homo<-plots1[which(plots1$chra==plots1$chrb),]
pos=with(homo,cbind(starta,startb))
#plot(pos, pch='.', col="#77777777")

chrlen = max(pos)
gridsize = ceiling(chrlen/4e4)
bandwidth = 4e4
den = bkde2D(pos, bandwidth=c(1,1)*bandwidth, gridsize=c(1,1)*gridsize)
den$fhat <- den$fhat + t(den$fhat)

with(den, image(x=x1, y=x2, z=log(fhat), col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE))

r <- as(den$fhat, "RasterLayer")
rmean <- focal(r, w=matrix(1/49, ncol=7, nrow=7), fun=mean)
meand<-diag(as(rmean,"matrix"))

plot(log(meand),log(homomean))



m = matrix(0, nrow=gridsize, ncol=gridsize)

for(i in 1:(gridsize-1)) {
  band = (row(m)==col(m)+i)
  print(band)
  m[band] = mean(den$fhat[band])
}

m = m + t(m)
diag(m) = mean(diag(den$fhat))
image(x=den$x1, y=den$x2, z=m^0.3, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE)


genomicDistance = den$x1 - min(den$x1)
averageInteractions = m[1,]
plot(genomicDistance, sqrt(averageInteractions), type ="l")

fhatNorm <- den$fhat/m
image(x=den$x1, y=den$x2, z=fhatNorm, col=colorRampPalette(c("white","blue"))(256), useRaster=TRUE)

cm <- cor(fhatNorm)
image(x=den$x1, y=den$x2, z=cm,zlim=c(-1,1), col=colorRampPalette(c("red", "white","blue"))(256), useRaster=TRUE)

library(Sushi)

#Sushi_data = data(package ='Sushi')
#data(list = Sushi_data$results[,3])
#Sushi_data$results[,3]

# makepdf = TRUE
# pdfname="figure1_label.pdf"
# 
# if( makepdf == TRUE)
# {
#   pdf ( pdfname ,  height =10, width=12)
#   
# }

chrom = "chr1"
chromstart       = 40000000
chromend         = 80000000

# layout (matrix(c(1,1,
#                  2,2,
#                  3,3,
#                  4,4,
#                  5,11,
#                  6,12,
#                  7,13,
#                  8,14,
#                  9,15,
#                  10,16),10,2,byrow=TRUE))

par(mar=c(0.4,2,0.3,3))
dog1<-read.table("log2.ratio.preL_2N.mm10_rep_timing.bedgraph")
dog1$V4<-dog1$V4-min(dog1$V4)
plotBedgraph(dog1, chrom, chromstart, chromend,color="#FF6666")
axis(side=2,las=2,tcl=.2)
mtext("Replication timing",side=2,line=2.5,cex=1,font=2)
legend("topleft",legend=c(""),text.font=2,bty='n')


