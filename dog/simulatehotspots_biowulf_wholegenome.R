
## Author: Qi Yu
## Contact: dryuqi@gmail.com

simulatehotspots_biowulf <- function (
  nI             = 5000, #number of cells
  #plotAnimations = FALSE,
  genome         = 'canFam3',
  #chrs           = 'chr38',
  #win            = 1000, #simulation a window size
  sam            = 1, #the number of random sample 
  dsb            = 5,
  tel_weight     = 0,
  tss_weight     = 0,
  h3k4_weight     = 0,
  target_hotspot  = "/data/yuq3/simulation/dog1_noC_merge_longest_peak_noRC"
  #useStrength    = TRUE,
  #useTel_dis     = TRUE,
  #repTimingData  = NULL,
  #IODtime        = 50
  )
  
{
  
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library("R.utils")
library(foreach)
#library(doParallel)
#registerDoParallel(cores = 24)
#library(MASS)

#args<-(commandArgs(TRUE))

doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -
                                                 min(x, na.rm=TRUE))}

#if(!file.exists(paste0('genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))){
  
#system(paste('grep -P ',chrs,'\t', ' genome.fa_test_CG_all.distribute.flip.bed >genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))
#system(paste('grep -P ',chrs,'\t', ' dog_tel_new.bed_CG_all_distance.txt.flip.bed >dog_tel_new.bed_CG_all_distance.txt.flip.bed_',chrs,sep=''))
#system(paste('grep -P ',chrs,'\t', ' tss_center_canfam3_nochrX.bed.flip.bed_CG_all_distance.txt.flip.bed >tss_center_canfam3_nochrX.bed.flip.bed_CG_all_distance.txt.flip.bed_',chrs,sep=''))
#}

#if(!file.exists(paste0('lz_h3k4_genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))){
#system(paste('grep -P ',chrs,'\t', ' lz_h3k4_genome.fa_test_CG_all.distribute.flip.bed >lz_h3k4_genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))
#}

#if(!file.exists(paste0('genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))){
#system(paste('grep -P ',chrs,'\t', 'genome.fa_test_CG_all.distribute.flip.bed >genome.fa_test_CG_all.distribute.flip.bed_',chrs,sep=''))
#}
file<-read.table('CGI-Cfamiliaris.bed.flip.bed',header=FALSE)
#file<-read.table(paste("genome.fa_test_CG_all.distribute.flip.bed_",chrs,sep=''),header=FALSE)
#file<-read.table(paste("tss_2k_noCGI.bed_",chrs,sep=''),header=FALSE)
#file<-subset(file,(file$V1==chrs))
names(file)<-c("chr","start","end","len","CpG","CG","CpG_content","GC_content")
file$cen<-as.integer(file$start+(file$end-file$start)/2)
#file$weight<-nchar(as.character(file$CG))
file$weight<-file$CpG
###############################################weight by tel, tss or h3k4######################################################

weit=weit1=weit2=weit3=""
new<-file

#if (tel_weight == TRUE){

#tel<-read.table(paste('dog_tel_new.bed_CG_all_distance.txt.flip.bed_',chrs,sep=''),header=FALSE)
#tel<-read.table('/data/yuq3/simulation/dog_tel_new.bed_CG_all_distance.txt.flip.bed',header=FALSE)
#tel<-subset(tel,tel$V1==chrs)
#names(tel)<-c("chr","start","end","dtel")
#tel$dtel<-doit(-(log(tel$dtel+1)))
#new<-join(file,tel,by=c("chr","start","end"))
#new$tel_weight2<-new$dtel*tel_weight
weit<-paste("_By_telD",tel_weight,sep="")
#}

#if (tel_weight == FALSE){
#  new$tel_weight<-0
#}

#if (tss_weight == TRUE){
#tss<-read.table('/data/yuq3/simulation/tss_center_canfam3_nochrX.bed.flip.bed_CG_all_distance.txt.flip.bed',header=FALSE)
#tss<-subset(tss,tss$V1==chrs)
#names(tss)<-c("chr","start","end","dtss")
#tss$dtss<-doit(-(log(tss$dtss+1)))
#new<-join(new,tss,by=c("chr","start","end"))
#new$tss_weight2<-new$dtss*tss_weight
weit1<-paste("_By_tssD",tss_weight,sep="")
#}

#if (tss_weight == FALSE){

#  new$tss_weight<-0
#}

#if (h3k4_weight==TRUE){
#h3k4<-read.table('/data/yuq3/simulation/dog1_h3k4me3_middle.bed_CG_all_distance.txt.flip.bed',header=FALSE)
#h3k4<-subset(h3k4,h3k4$V1==chrs)
#names(h3k4)<-c("chr","start","end","dh3k4")
#h3k4$dh3k4<-doit(-(log(h3k4$dh3k4+1)))
#h3k4$h3k4_weight[which(h3k4$dh3k4<3000)]<-1
#h3k4$h3k4_weight[which(h3k4$dh3k4>=3000)]<-0.01
#h3k4 <- h3k4[ which(h3k4$h3k4_weight==1),]
#h3k4$dh3k4<-doit(-(log(h3k4$dh3k4+1)))
#new<-join(new,h3k4,by=c("chr","start","end"))
#new$h3k4_weight2<-new$dh3k4*h3k4_weight
weit2<-paste("_By_h3k4D",h3k4_weight,sep="")
#}

#if (h3k4_weight == FALSE){

#  new$h3k4_weight<-0
#}

#if (  tel_weight == FALSE && tss_weight == FALSE && h3k4_weight == FALSE){

#  new$newweight<-new$weight
#  weit3<-"noneweight"
#}

weinew<-paste(weit,weit1,weit2,weit3,sep="")
new$newweight<-new$weight
######################################################################################################

if (genome == 'canFam3'){
  genomeSizes         <- read.table('/data/yuq3/genome/canFam3/genome.fa.fai',
                                    header=FALSE)
  names(genomeSizes)  <- c('cs','size','cum','xa','xb')
  genomeSizes<-subset(genomeSizes,genomeSizes$cs != "chrX" & genomeSizes$cs != "chrM")

  #chromSize           <- genomeSizes$size[genomeSizes$cs == chrs]
}

#dsb<-round(chromSize/4000000)
#################################################################################################
#file<-new
#nRow = chromSize
options(scipen=10000)
if(!file.exists(paste(nI,"simulation_CGI_binorm_750",weinew,'read',sam,'dsb.bed',sep="_"))){
 
#for (i in 1:nI){
for (ch in genomeSizes$cs){  
  chromSize=genomeSizes$size[genomeSizes$cs == ch]
  dsb<-round(chromSize/10000000)
  file<-subset(new,new$chr==ch)
  for (i in 1:nI) {
  randPositions <- sample((file$cen),dsb,prob=file$newweight,replace=TRUE)
  for (x in randPositions){

    rpos<-as.integer(rnorm(5000,mean=-750,sd=400))
    rpos2<-as.integer(rnorm(5000,mean=750,sd=400))
    randReads<- sample(rpos,sam)
    randReads2<- sample(rpos2,sam)
    for (y in c(randReads,randReads2)){
      if (y+50<0){
        write(paste(ch,x+y,x+y+50,"reverse","+",sep="\t"),file=paste(nI,"simulation_CGI_binorm_750",weinew,'read',sam,'dsb.bed',sep="_"),append=TRUE)
      }
      if (y>0){
        write(paste(ch,x+y,x+y+50,"forward","-",sep="\t"),file=paste(nI,"simulation_CGI_binorm_750",weinew,'read',sam,'dsb.bed',sep="_"),append=TRUE)
      }
      
    }
  }
 } 
}
}
###############################################count overlap and calculate correlation##########################################
readf<-function(x)
{
  a<-read.table(x, header=TRUE)

  #names(a)=c("chr","start","end","tags")
  #a <- a[which(a$cs==chrs),]
  a<-a[,c("cs","from","to","strength")]
  names(a)=c("chr","start","end","tags")
  
  return(a)
}

##system(paste('grep ',chrs,' 120229_0179_1_NoIndex_FS_021712_dog1_dmc1.canfam3_flip.ssDNA_type1.bed.60_60.bed > ',chrs,'_allreads.bed',sep=''))
##sumr<-countLines(paste(chrs,'_allreads.bed',sep=''))
ot<-readf(target_hotspot)

options(scipen=10000)

system(paste('sh bed_biowulf_wholegenome.sh ',nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb.bed ',' ',target_hotspot,sep=''))

##system(paste('perl ./calcStrengthAndRecenterHotspots.pl --hs dog_peak_',chrs,' --frag ',nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed',' --noRC ',' --out ',nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed_peak.bed',sep=''))
#suml<-countLines(paste(nI,'_simulation_CG_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed',sep=''))
#system(paste('coverageBed -a ',nI,'_simulation_CG_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed', ' -b dog_peak_',chrs,' > ',nI,'_simulation_CG_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed_peak.bdg',sep=''))
  
##ol<-readf(paste(nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb_',chrs,'_.bed.bed_peak.tab',sep=''))
##mergeo<-data.frame(join_all(list(ol,ot),by=c("chr","start","end")))
##names(mergeo)=c("chr","start","end","simu","real")
  
##b<-cor(mergeo$simu,mergeo$real,method="spearman")
  
TP<-countLines(paste(nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb.bed_overlap.bed',sep=''))
FP<-countLines(paste(nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb.bed_nonoverlap.bed',sep=''))
  
##b<-cor(mergeo$simu,mergeo$real,method="spearman")
#write(paste(nI,sam,dsb,chrs,tel_weight,tss_weight,h3k4_weight,TP,FP,sep="\t"),file="grid_search_TP_FP_result.txt",append=TRUE)
  

################################read coverage###########################3

readc<-function(x)
{
  a<-read.table(x, header=FALSE)
  names(a)<-c("chr","start","end","cover")
  a$cover <- a$cover*1000000
  a<-a[,c("chr","start","end","cover")]
  return(a)
}

dog1<-readc("/data/yuq3/simulation/grid_search_whole_genome/CpG/120229_0179_1_NoIndex_FS_021712_dog1_dmc1.canfam3_flip.ssDNA_type1.bed.60_60.bed_100kb_4line.bdg")

#dog1<-subset(dog1,dog1$chr==chrs)
################read domain###############################################################
#dog1d<- read.table("dog_domain_manual_3types_10M_nopeak_dog1_sort.bdg")
#names(dog1d)=c("chr","start","end","loc","Domain","subtract")

#dog1d$start<-dog1d$start/1000000
#dog1d$end<-dog1d$end/1000000
#dogd<-subset(dog1d,dog1d$chr==chrs)

#options(scipen=10000)
################################read 100kb coverage############################################
sim<- readc(paste(nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb.bed_100kb.bedgraph',sep=''))

#merge<-data.frame(join_all(list(dog1,dog2,dog3),by="ID"))

#names(merge)=c("Dog1","position","Dog2","Dog3")
#new<-melt(merge,id=c("position"))
#doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -
               ##                                  min(x, na.rm=TRUE))}

final<-data.frame(join_all(list(dog1,sim),by=c("chr","start","end")))
names(final)=c("chr","start","end","real","simu")
final$coverage=doit(final$real)
final$simulation=doit(final$simu)
fm<-melt(final[,c("start","coverage","simulation")], id=c("start"))
#############################################################################3
getRMSE = function(m, o){
  NAok <- !is.na(m*o)
  retRMSE <- sqrt(mean((m[NAok] - o[NAok])^2))
  return(retRMSE)
}

######################################################################
### Sum square Error Function
getSSE = function(m, o){
  NAok   <- !is.na(m*o)
  retSSE <- sum((m[NAok] - o[NAok])^2)
  return(retSSE)
}

R2   <- cor(final$coverage,final$simulation,use='complete.obs')
RMSE <- getRMSE(final$coverage,final$simulation)
SSE  <- getSSE(final$coverage,final$simulation)

write(paste(nI,sam,dsb,tel_weight,tss_weight,h3k4_weight,TP,FP,R2,RMSE,SSE,sep="\t"),file="grid_search_TP_FP_result.txt",append=TRUE)


##write(paste(nI,sam,dsb,chrs,weinew,b,TP,FP,R2,RMSE,SSE,sep="\t"),file="large_scale_correlation_result_new.txt",append=TRUE)
#############################################################################3

##d=ggplot()+
##  geom_line(data=fm,aes(x=(fm$start)/1000000,y=fm$value,color=variable),size=1)+
##  geom_rect(data=dogd,aes(xmin=dogd$start,xmax=dogd$end,ymin=-0.05,ymax=0),fill="gray",color="gray")+
  #geom_line(data=dog1,aes(x=dog1$ID,y=doit(dog1$V1)),size=1)+
##  theme_bw()+
##  ylab("Normalized coverage")+
  #scale_y_log10()+
  #scale_x_continuous(limits=c(0,38), breaks =c(2,20,38), labels = c("-2k","Hotspots","2k"))+
  #annotate("text", x = 0.8, y = 5, label = c("human"),colour="blue"))+
##  ggtitle(label = sprintf(" %s cells; DSB = %i ; read = %s; weighted=%s;\n R2 = %4.2f; RMSE = %5.4f; SSE = %5.4f",
##                          nI,
##                          dsb,
##                          sam,
##                          weinew,
##                          R2,
##                          RMSE,
##                          SSE)) + 
  
##  guides(color=guide_legend( keywidth=3, keyheight=3, title=""))+
  #theme(legend.position=c(0.5,0.6))+
  #scale_color_manual(values=c("#cc0000","#009900"),breaks=c("Dog1","Simulated"))+
##  theme(text = element_text(size=30,family="sans"),strip.background=element_rect(fill="white"))+
##  xlab(" ")
##ggsave(plot=d,filename=paste(nI,'_simulation_CGI_binorm_750_',weinew,'_read_',sam,'_dsb_',chrs,'.bed_100kb','_coverabe_compare.tiff',sep=""),height=6,width=16)

return ()
}

