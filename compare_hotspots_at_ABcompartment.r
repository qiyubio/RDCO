## Author: Qi Yu
## Contact: dryuqi@gmail.com

readf<-function(x,y)
{
  a<-read.table(x, header=FALSE)
  names(a)=c("chr","start","end","strength","a","b","c","EI")

  a$type<-y 
  return(a)
}

readfa<-function(x,y)
{
  a<-read.table(x, header=FALSE)
  names(a)=c("chr","start","end","a","b","c","EI")
  
  a$type<-y 
  return(a)
}


data1<- readf("B6_hotspots_eigenvector.bed","Hotspots")
data2<- readfa("random_hotspots_noX_noY_lep_EI.bed","random")

data2$strength<-0

first<-rbind(data1,data2)


d<-first[which(first$type=="Hotspots"),]
n<-first[which(first$type=="random"),]
wilcox.test(d$EI,n$EI,alternative="greater")


first %>%
  group_by(name,type) %>%
  summarize(median(EI))

first$strength_nor<-first$strength/0.3
firstc<-c("name","type","strength_nor")
DMC1s<-first[firstc]

new1<-rbind(DMC1s,Rbs)

d1=ggplot(first,aes(x=type,y=EI,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.1, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  #scale_x_discrete(expand=c(0.2,0))+
  annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  
  ylab("EI")+
  xlab("")+
  #scale_y_log10()+
  #ylim(0,3)+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  #scale_colour_discrete(name = "")+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  #scale_color_manual(values=c("#cc0000","#009900","#ff9900"),breaks=c("Dog1","Dog2","Dog3"))+
  ggtitle("")


ggsave(plot=d1,filename=paste("EI_hotspots_random_Lep.tiff",sep=""),height=5,width=4,dpi=200)
ggsave(plot=d1,filename=paste("Extra peaks_random_peak.pdf",sep=""),height=5,width=4,dpi=200)

################################################################################################
data1<- readf("30000_simulation_CG_logweight_binorm_750_noneweight_read_1_dsb_8_autosome_.bed_macs2_peaks.xls.bed_nonoverlap_peak_at_nondomain_short.bed.coverage.bed_FPKM.bed","Extra peaks")
data2<- readf("random_nondomain_9500.bed.coverage.bed_FPKM.bed","random")

first<-rbind(data1,data2)
first$name<-"DMC1"

table<-read.table("mark4_cleaned_canFam3.1.RR_CM.bed.flip.bed_30000_simulation_CG_logweight_binorm_750_noneweight_read_1_dsb_8_autosome_.bed_macs2_peaks.xls.bed_nonoverlap_peak_at_nondomain_short.bed_only_per_chr.bed",sep="\t")
names(table)=c("chr","start","end","Rb")
table$type<-"Extra peaks"


table2<-read.table("mark4_cleaned_canFam3.1.RR_CM.bed.flip.bed_random_nondomain_9500.bed_only_per_chr.bed",sep="\t")
names(table2)=c("chr","start","end","Rb")
table2$type<-"random"

d<-first[which(first$type=="Extra peaks"),]
n<-first[which(first$type=="random"),]
wilcox.test(d$strength,n$strength,alternative="greater")

d<-second[which(second$type=="Extra peaks"),]
n<-second[which(second$type=="random"),]
wilcox.test(d$strength,n$strength,alternative="greater")

second<-rbind(table,table2)
second$name<-"RR"

second %>%
  group_by(name,type) %>%
  summarize(median(Rb))

second$strength_nor<-second$Rb/0.897
selectc<-c("name","type","strength_nor")
Rbs<-second[selectc]
##
first %>%
  group_by(name,type) %>%
  summarize(median(strength))

first$strength_nor<-first$strength/0.29
firstc<-c("name","type","strength_nor")
DMC1s<-first[firstc]


new2<-rbind(DMC1s,Rbs)

d2=ggplot(new2,aes(x=name,y=strength_nor,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  #scale_x_discrete(expand=c(0.2,0))+
  annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  ggtitle("nonDomain")+
  ylab("Relative strength")+
  xlab("nonDomain")+
  scale_y_log10()+
  #ylim(0,3)+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  #scale_colour_discrete(name = "")+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  #scale_color_manual(values=c("#cc0000","#009900","#ff9900"),breaks=c("Dog1","Dog2","Dog3"))+
  ggtitle("")


ggsave(plot=d2,filename=paste("Extra peaks_random_peak_nonDomain.tiff",sep=""),height=5,width=4,dpi=200)
ggsave(plot=d2,filename=paste("Extra peaks_random_peak_nonDomain.pdf",sep=""),height=5,width=4,dpi=200)

#########################################################################################

data1<- readf("30000_simulation_CG_logweight_binorm_750_noneweight_read_1_dsb_8_autosome_.bed_macs2_peaks.xls.bed_nonoverlap_peak_at_domain.bed.short.bed.coverage.bed_FPKM.bed","Extra peaks")
data2<- readf("random_domain_5000.bed.coverage.bed_FPKM.bed","random")

first<-rbind(data1,data2)
first$name<-"DMC1"

table<-read.table("mark4_cleaned_canFam3.1.RR_CM.bed.flip.bed_30000_simulation_CG_logweight_binorm_750_noneweight_read_1_dsb_8_autosome_.bed_macs2_peaks.xls.bed_nonoverlap_peak_at_domain.bed.short.bed_only_per_chr.bed",sep="\t")
names(table)=c("chr","start","end","Rb")
table$type<-"Extra peaks"


table2<-read.table("mark4_cleaned_canFam3.1.RR_CM.bed.flip.bed_random_domain_5000.bed_only_per_chr.bed",sep="\t")
names(table2)=c("chr","start","end","Rb")
table2$type<-"random"

d<-first[which(first$type=="Extra peaks"),]
n<-first[which(first$type=="random"),]
wilcox.test(d$strength,n$strength,alternative="greater")

d<-second[which(second$type=="Extra peaks"),]
n<-second[which(second$type=="random"),]
wilcox.test(d$strength,n$strength,alternative="greater")

second<-rbind(table,table2)
second$name<-"RR"

second %>%
  group_by(name,type) %>%
  summarize(median(Rb))

second$strength_nor<-second$Rb/1.26
selectc<-c("name","type","strength_nor")
Rbs<-second[selectc]
##
first %>%
  group_by(name,type) %>%
  summarize(median(strength))

first$strength_nor<-first$strength/0.59
firstc<-c("name","type","strength_nor")
DMC1s<-first[firstc]


new3<-rbind(DMC1s,Rbs)

d3=ggplot(new3,aes(x=name,y=strength_nor,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  #scale_x_discrete(expand=c(0.2,0))+
  annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  
  ylab("Relative strength")+
  xlab("Domain")+
  scale_y_log10()+
  #ylim(0,3)+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  #scale_colour_discrete(name = "")+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  #scale_color_manual(values=c("#cc0000","#009900","#ff9900"),breaks=c("Dog1","Dog2","Dog3"))+
  ggtitle("")


ggsave(plot=d3,filename=paste("Extra peaks_random_peak_Domain.tiff",sep=""),height=5,width=4,dpi=200)
ggsave(plot=d3,filename=paste("Extra peaks_random_peak_Domain.pdf",sep=""),height=5,width=4,dpi=200)


yl<-max(new1$strength_nor,new2$strength_nor,new3$strength_nor)
ym<-min(new1$strength_nor,new2$strength_nor,new3$strength_nor)+0.01

d1=ggplot(new1,aes(x=name,y=strength_nor,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  #scale_x_discrete(expand=c(0.2,0))+
  annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  ylab("Relative strength")+
  xlab("All")+
  scale_y_log10(limits=c(ym,yl))+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  ggtitle("")

d2=ggplot(new2,aes(x=name,y=strength_nor,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  xlab("nonDomain")+
  scale_y_log10(limits=c(ym,yl))+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="top",)+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  ggtitle("")
  
d3=ggplot(new3,aes(x=name,y=strength_nor,colour=type),)+
    theme_bw()+
    geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
    geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
    
    xlab("Domain")+
    scale_y_log10(limits=c(ym,yl))+
    theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),axis.title.y=element_blank(),axis.text.y=element_blank(),legend.position="top",)+
    guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
    scale_color_manual(values=c("#C35817","grey"),name="")+
    ggtitle("")
  
tiff(filename =paste("Extra peaks_random_peak_compare.tiff"), width = 22, height = 10,units="cm", res = 300)
grid.arrange(d1,d2,d3, nrow=1,ncol=3,widths = c(1,0.8,0.8))
dev.off() 

pdf("Extra peaks_random_peak_compare.pdf", width = 8, height = 4)
grid.arrange(d1,d2,d3, nrow=1,ncol=3,widths = c(1,0.8,0.8))
dev.off() 

#######################################################################################
readf<-function(x,y)
{
  a<-read.table(x, header=FALSE)
  names(a)=c("chr","start","end","name","strength")
  
  a$type<-y 
  return(a)
}

data1<- readf("30000_simulation_CG_mm10_binorm_750_noneweight_read_merge_peak.bed_nonoverlap_hotspots.bed.FPKM.bed","Extra peaks")
data2<- readf("random_simulation_peak_nonoverlap_hotspots.bed.FPKM.bed","random")

first<-rbind(data1,data2)
first$name<-"DMC1"



d<-first[which(first$type=="Extra peaks"),]
n<-first[which(first$type=="random"),]
wilcox.test(d$strength,n$strength,alternative="greater")


first %>%
  group_by(name,type) %>%
  summarize(median(strength))

first$strength_nor<-first$strength/0.223
firstc<-c("name","type","strength_nor")
DMC1s<-first[firstc]


new2<-DMC1s

d2=ggplot(new2,aes(x=name,y=strength_nor,colour=type),)+
  theme_bw()+
  geom_violin(scale = "width",size=1.5,trim = TRUE,position=position_dodge(0.85))+
  geom_boxplot(width=0.4, notch=TRUE, position=position_dodge(0.85),size=1.5)+
  #scale_x_discrete(expand=c(0.2,0))+
  annotation_logticks(sides="l")+
  #geom_dotplot(binaxis = "y", stackdir = "center",alpha=0.5)+
  ggtitle("nonDomain")+
  ylab("Relative strength")+
  xlab("nonDomain")+
  scale_y_log10()+
  #ylim(0,3)+
  theme(text = element_text(size=20,family="sans"),axis.text.x=element_text(size=20,color="black"),legend.position="top",)+
  #scale_colour_discrete(name = "")+
  guides(colour=guide_legend( keywidth=2, keyheight=2, title=""))+
  scale_color_manual(values=c("#C35817","grey"),name="")+
  #scale_color_manual(values=c("#cc0000","#009900","#ff9900"),breaks=c("Dog1","Dog2","Dog3"))+
  ggtitle("")


ggsave(plot=d2,filename=paste("prdm9ko_Extra_peaks_random_peak_nonDomain.tiff",sep=""),height=5,width=4,dpi=200)
ggsave(plot=d2,filename=paste("prdm9ko_Extra_peaks_random_peak_nonDomain.pdf",sep=""),height=5,width=4,dpi=200)

