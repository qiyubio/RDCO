## Author: Qi Yu
## Contact: dryuqi@gmail.com

library(ggplot2)
library(wesanderson)
args<-(commandArgs(TRUE))
library(tidyr)
library(tibble)
library(reshape2)

#filename = args[1]
#title = args[2]

readfile<-function(filename,title=NA){
  homo<-read.table(filename)
  homo<-homo[c(2,5,7)]
  names(homo)<-c("start","count","end")
  homo$distance<-abs(homo$end-homo$start)/1000
  
  d <- homo$distance
  cs <- sum(homo$count)
  #100k to 2000k
  intvals <- 1:40
  props <- sapply(intvals, function(i){
    sum(homo$count[(d <= i * 100) & (d > (i-1) * 100)])/cs
  })
  
  v <- as.character(c(1,100*intvals)); lv <- length(v)
  colnm <- paste0(c(paste0(v[-1*lv], "-",v[-1])[-lv]), "kb")
  names(props) <- colnm
  propdf <- melt(props)
  propdf$type<-title
  propdf <- rownames_to_column(propdf, var="Bin")
  return(propdf)
}#end of function readfile

table<-readfile("/home/qiyu/data/hiC/hichipper_lep_analysis/190110_9018_1_CGTACG_GC_121718_HiChip_B6_4C_H3K4me3.intra.loop_counts.bedpe_merge_at_B6_wt_DMC1_peaks_no_blacklist.bed","4C_1HT")
table2<-readfile("/home/qiyu/data/hiC/hichipper_lep_analysis/190110_9018_1_CGATGT_GC_121718_HiChip_B6_Lep_H3K4me3.intra.loop_counts.bedpe_merge_at_B6_wt_DMC1_peaks_no_blacklist.bed","Lep_1HT")

# table<-table[which(table$distance<2000),]
# table2<-table2[which(table2$distance<2000),]
# 
# new<-table(cut(table$distance, seq(0, 2000, by=100), include.lowest=TRUE))
# new2<-table(cut(table2$distance, seq(0, 2000, by=100), include.lowest=TRUE))
# 
# tablea<-split(table, cut(table$distance, seq(0, 2000, by=100), include.lowest=TRUE))
# tableb<-split(table2, cut(table2$distance, seq(0, 2000, by=100), include.lowest=TRUE))
# 


# Dummy calls to remove note
Sample <- ""
Bin <- ""
Value <- ""

intvals <- 1:40
v <- as.character(c(1,100*intvals)); lv <- length(v)

propdf <- rbind(table,table2)
colnm <- paste0(c(paste0(v[-1*lv], "-",v[-1])[-lv]), "kb")
propdf$Bin<-factor(propdf$Bin,levels=colnm,ordered=TRUE)
propdf<-propdf[order(propdf$Bin),]

p1 <- ggplot(propdf, aes(x=Bin, y=value, group=type, color=type)) + geom_line(size=2) +
  theme_bw()+
  theme(axis.text.x = element_text(hjust=1, angle = 45), 
        text = element_text(size=25),
        legend.title=element_blank()) +
  labs(title = "Proportional PET counts per distance",
       x = "Binned Loop Distance", y = "Proportional Pet Counts") +
  annotation_logticks(sides="l")+
  scale_y_log10()



ggsave(plot=p1,filename="intra_noself_1HT_4Mb_4C_Lep.png",height=7,width=14,dpi=200)

#new2$type<-factor(new2$type,levels=c("PreL_Boris","P_Boris","P_CTCF","Z_Boris","Z_CTCF","P_Boris_DSB","P_CTCF_DSB"))
options(scipen=10000)

new4<-new2[which(new2$distance<1000),]
d2<-ggplot(new4,aes(x=distance,fill=type))+
  #geom_density(linesize=2)+
  geom_histogram(color="black")+
  facet_wrap( ~ type, ncol=1)+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("loop Length (kb)")+
  #scale_y_log10()+
  annotation_logticks(sides="b")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="1HT_intra_noself_4C_Lep.png",height=7,width=10,dpi=200)

d2<-ggplot(new2,aes(x=count,color=type))+
  #geom_density(linesize=2)+
  stat_ecdf(size=2)+
  facet_wrap( ~ type, ncol=1)+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("PET")+
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks(sides="lb")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="PET_count_noself_4C_Lep.png",height=7,width=10,dpi=200)

new3<-new4[which(new4$count>4),]
d2<-ggplot(new3,aes(x=distance,fill=type))+
  #geom_density(linesize=2)+
  geom_histogram(color="black")+
  facet_wrap( ~ type, ncol=1)+
  ggtitle(">=5")+
  #facet_grid(type ~ .)+
  #geom_freqpoly()+
  #geom_density(size=1.5)+
  #scale_fill_manual(values=c(rev(wes_palette("Rushmore",5)),"blue","purple"),name="")+
  theme_bw()+
  xlab("loop Length (kb)")+
  #scale_x_log10()+
  annotation_logticks(sides="b")+
  theme(text = element_text(size=25),legend.title=element_blank(),legend.justification=c(1,01))
#scale_x_log10()+scale_y_log10()+
#scale_x_log10()+scale_y_log10()+

ggsave(plot=d2,filename="1HT_intra_noself_4C_Lep_5_a.png",height=7,width=10,dpi=200)


