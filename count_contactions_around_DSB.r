## Author: Qi Yu
## Contact: dryuqi@gmail.com


library(ggplot2)
library(HiTC)
library(scales)
library(dplyr)
library(plyr)
library(tidyr)


#end of function plotIntraDistnew

readfile<-function(filename,title=NA,genome="mm10"){
  homo<-read.table(filename)
  names(homo)<-c("hchr","hstart","hend","chr1","start1","end1","chr2","start2","end2")
  homo$distance<-abs(homo$end2-homo$end1)
  homo$type<-title
  homo<-homo[,c("hchr","distance","type")]
  return(homo)
  #plotIntraDistnew(homo,winsize=40000,title=title,genome=genome)
}#end of function readfile

#readfile("/home/yuqi/Downloads/GSM2745898_17MAY11_HSC02RPAB-PE50_CCHiC-HeLa-NS-R1_23-05-2011_kittlere.5_hg19.validPair.txt.gz_count.txt",title="hela_test",genome="hg19")
#readfile("/mnt/RDCO/shared/Gang/sperm/merged_nodups.txt_count.txt",title="sperm",genome="mm10")
lep<-readfile("random_8k_noXY.bed_spermatid.txt",title="random",genome="mm10")
prel<-readfile("B6_hotspots_8k.bed_spermatid.txt",title="DSB",genome="mm10")

merge<-rbind(lep,prel)

d=ggplot(merge)+
  geom_density(aes(x=abs(distance),color=type),size=1)+
  xlab("Genomic Distance")+
  ylab("Density")+
  theme_bw()+
  theme(text = element_text(size=25))+
  #ggtitle(title)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="b")
ggsave(plot=d,filename=paste("read_pair_distance_from_DSB_or_not_spermatid.png", sep=''),height=10,width=16)

merge$distance<-abs(merge$distance)
names(merge)<-c("chr","distance","type")
mergea<-merge[which(merge$distance >=1000 & merge$distance <2000000),]
mergeb<-merge[which(merge$distance >=2000000),]

aggdata <-with(mergea,table(chr,type))
aggdatab <-with(mergeb,table(chr,type))

aggdata1<-data.frame(aggdata)
aggdatab1<-data.frame(aggdatab)

aggdata1<-aggdata1[which(aggdata1$chr!="chrM" & aggdata1$chr!="chrY"),]
aggdatab1<-aggdatab1[which(aggdatab1$chr!="chrM" & aggdata1$chr!="chrY"),]

merge2M<-join(aggdata1,aggdatab1,by=c("chr","type"))
names(merge2M)<-c("chr","type","less2M","more2M")
merge2M$ratio<-merge2M$more2M/merge2M$less2M

#merge2M$type<-factor(merge2M$type,levels=c("preLep","Lep","zygotene","pachytene","diplotene","spermatid"))

d=ggplot(merge2M,aes(x=type,y=ratio))+
  geom_violin(size=1.5)+
  geom_boxplot(width=0.2, notch=TRUE, position=position_dodge(0.9),size=1.5)+
  ylab("Ratio of interactions (>=2M) vs (<2M)")+
  xlab(" ")+
  theme_bw()+
  theme(text = element_text(size=20),axis.text.x = element_text(angle = 50, hjust = 1))
#ggtitle(title)+
#scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
#scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
#annotation_logticks(sides="b")
ggsave(plot=d,filename=paste("ratio_around_DSB_or_random.png", sep=''),height=6,width=6)


mergea20k<-merge[which(merge$distance >=1000 & merge$distance <20000),]
mergea40k<-merge[which(merge$distance >=20000 & merge$distance <40000),]
mergea60k<-merge[which(merge$distance >=40000 & merge$distance <60000),]
mergea80k<-merge[which(merge$distance >=60000 & merge$distance <80000),]
mergea100k<-merge[which(merge$distance >=80000 & merge$distance <100000),]
mergea120k<-merge[which(merge$distance >=100000 & merge$distance <120000),]
mergea140k<-merge[which(merge$distance >=120000 & merge$distance <140000),]
mergea160k<-merge[which(merge$distance >=140000 & merge$distance <160000),]


aggdata20 <-with(mergea20k,table(chr,type))
aggdata40 <-with(mergea40k,table(chr,type))
aggdata60 <-with(mergea60k,table(chr,type))
aggdata80 <-with(mergea80k,table(chr,type))
aggdata100 <-with(mergea100k,table(chr,type))
aggdata120 <-with(mergea120k,table(chr,type))
aggdata140 <-with(mergea140k,table(chr,type))
aggdata160 <-with(mergea160k,table(chr,type))

aggdata20d<-data.frame(aggdata20)
aggdata40d<-data.frame(aggdata40)
aggdata60d<-data.frame(aggdata60)
aggdata80d<-data.frame(aggdata80)
aggdata100d<-data.frame(aggdata100)
aggdata120d<-data.frame(aggdata120)
aggdata140d<-data.frame(aggdata140)
aggdata160d<-data.frame(aggdata160)


removeM_Y<-function(test){
  test<-test[which(test$chr!="chrM" & test$chr!="chrY" & test$chr!="chrX"),]
  return(test)
}

aggdata20d<-removeM_Y(aggdata20d)
aggdata40d<-removeM_Y(aggdata40d)
aggdata60d<-removeM_Y(aggdata60d)
aggdata80d<-removeM_Y(aggdata80d)
aggdata100d<-removeM_Y(aggdata100d)
aggdata120d<-removeM_Y(aggdata120d)
aggdata140d<-removeM_Y(aggdata140d)
aggdata160d<-removeM_Y(aggdata160d)


mergeallM<-join_all(list(aggdata1,aggdatab1,aggdata20d,aggdata40d,aggdata60d,aggdata80d,
                         aggdata100d,aggdata120d,aggdata140d,aggdata160d),by=c("chr","type"))

names(mergeallM)<-c("chr","type","less2M","more2M","20k","40k","60k","80k","100k","120k","140k","160k")

new<-sweep(mergeallM[,c(5,6,7,8,9,10,11,12)],mergeallM[,"less2M"]+mergeallM[,"more2M"],MARGIN=1,"/")

newmerge<-cbind(mergeallM[,c(1,2)],new)


cols.num <- c("20k","40k","60k","80k","100k","120k","140k","160k")

newmerge[cols.num] <- sapply(newmerge[cols.num],as.numeric)
sapply(newmerge, class)

newmerge %>%
  select(type,"20k","40k","60k","80k","100k","120k","140k","160k") %>% 
  #group_by(type) %>% 
  summarise_at(funs(wilcox.test(.[type=="DSB"],.[type=="random"])),.vars = newmerge$'20k')


wilcox.test(newmerge$'20k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'40k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'60k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'80k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'100k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'120k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'140k' ~ type, data = newmerge, alternative = "greater")$p.value
wilcox.test(newmerge$'160k' ~ type, data = newmerge, alternative = "greater")$p.value


newmerge<-melt(newmerge,id=c("chr","type"))
#newmerge$type<-factor(merge2M$type,levels=c("preLep","Lep","zygotene","pachytene","diplotene","spermatid"))

d=ggplot(newmerge,aes(x=variable,y=value,color=type))+
  #geom_violin(size=1.5)+
  geom_boxplot(width=0.8, notch=TRUE, position=position_dodge(0.9),size=0.5)+
  ylab("Propotion of short distance interactions")+
  xlab(" ")+
  scale_y_log10()+annotation_logticks(sides="l")+
  theme_bw()+
  theme(text = element_text(size=20),axis.text.x = element_text(angle = 50, hjust = 1))
#ggtitle(title)+
#scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
#scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
#annotation_logticks(sides="b")
ggsave(plot=d,filename=paste("ratio_around_each_stage_DSB_ornot2.png", sep=''),height=6,width=10)
