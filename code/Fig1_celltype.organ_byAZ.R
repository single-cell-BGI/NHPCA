####################   Rscript Fig1_celltype.organ_byAZ.R outdir infile date

library(RColorBrewer)
library(dplyr)
library(ggplot2)

args <- commandArgs(T)
outdir <- args[1]
infile <- args[2]
date <- as.character(args[3])


dir.create(outdir)
setwd(outdir)
info <- read.table(infile,header=T,sep="\t", stringsAsFactors=F,row.names=1)
info$organ <- gsub("_"," ",info$organ)


col_organ <- unique(c('#aec7e8','#98df8a','#c5b0d5', '#c49c94', '#f7b6d2',brewer.pal(9,"Pastel1")[c(1,2,4,5,7)],brewer.pal(12,"Set3")[c(1:8,10,11)],brewer.pal(8,"Set2")[c(3,5)],brewer.pal(12,"Paired")[c(1,2,5,7,9)],brewer.pal(8,"Accent")[c(1:4)],brewer.pal(11,"BrBG")[c(2,3,8:10)],brewer.pal(11,"PiYG")[c(2,3,9,10,11)],brewer.pal(11,"RdYlBu")[c(3,4,5,8,9)]))  ##45
col_organ <- colorRampPalette(col_organ)(length(table(info$organ)))[1:length(table(info$organ))]
set.seed(16)
col_organ <- sample(col_organ)

#### celltype
da <- data.frame(celltype=names(table(info$celltype)),value=c(paste0("0",1:9),10:length(table(info$celltype))))
rownames(da) <- da[,1]
info$celltype_lable <- da[as.character(info$celltype),"value"]

best_color<- unique(c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFE4E1","#0000A6","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#1B4400","#4FC601","#3B5DFF","#BA0900","#FF2F80","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#DDEFFF","#7B4F4B","#A1C299","#0AA6D8","#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF","#0067A5","#63FFAC","#F38400","#A1CAF1", "#C2B280","#848482","#E68FAC", "#F99379", "#604E97","#F6A600", "#B3446C","#DCD300","#882D17", "#8DB600","#654522", "#E25822", "#2B3D26", "#191970","#000080","#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493", "#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F", "#00FF00","#76EE00","#FFF68F","Yellow1","Gold1", "DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24", "#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347", "#FF4500","#FF6EB4","#EE30A7","#8B008B","green",brewer.pal(9,"Set1")[-9],brewer.pal(12,"Paired")[c(3,4,6,8,10,12)],brewer.pal(8, "Dark2"),'#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf' ,'#ffbb78','#ff9896')) ## unique 113
col_celltype <- colorRampPalette(best_color)(length(table(info$celltype)))[1:length(table(info$celltype))] ##113

set.seed(50)
col_celltype <- sample(col_celltype,length(col_celltype))
col_celltype <- col_celltype[order(as.numeric(levels(info$celltype_lable)))]
names(col_celltype) <- levels(info$celltype)


############################################ Fig1.1_Fig1A get umap_organ graph
da <- data.frame(organ=names(table(info$organ)),class=c(paste0("0",1:9),10:length(table(info$organ))),number=as.numeric(table(info$organ)))
write.table(cbind(da,col_organ),paste0("UMAP_organ_color.csv"),col.names=T,row.names=F,sep=",",quote=F)
rownames(da) <- da[,1]
info$organ_lable <- da[as.character(info$organ),"class"]
#organ lable
info %>%
  dplyr::group_by(organ_lable) %>%
  summarize(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2)) -> centers
#info$organ <- paste0(info$organ,".",info$organ_lable)

p1 <- ggplot(data=info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=0.1, size=0.00000001,aes(colour=organ))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_organ)+
	theme_classic()+
	labs(x='UMAP1',y='UMAP2',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank())+
	theme(legend.text=element_text(colour="black",size=12),legend.position="right",plot.title=element_text(hjust=0.5))+
  geom_point(data=centers, mapping=aes(x=UMAP_1, y=UMAP_2), size=0, alpha=0) +
  geom_text(data=centers, mapping=aes(label=organ_lable), size=4, fontface="plain")

pdf(paste0("UMAP_organ_",date,".pdf"),width=17,height=8)
print(p1)
dev.off()
png(paste0("UMAP_organ_",date,".png"),width=1100,height=500)
print(p1)
dev.off()


##################  umap_organ graph only legend
only.info <- info[!duplicated(info$organ),]
p1 <- ggplot(data=only.info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=1, size=6,aes(colour=organ))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_organ)+
	theme_classic()+
	labs(x='',y='',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),axis.line = element_blank())+
	theme(legend.text=element_text(colour="black",size=12),plot.title=element_text(hjust=0.5))+
  geom_point(data=centers, mapping=aes(x=UMAP_1, y=UMAP_2), size=0, alpha=0) +
  geom_text(data=centers, mapping=aes(label=organ_lable), size=3, fontface="plain")

p2 <- p1+theme(legend.position="right")
pdf(paste0("UMAP_organ_onlyLegend",date,".pdf"),width=15,height=8)
print(p2)
dev.off()
p2 <- p1+theme(legend.position="none")
pdf(paste0("UMAP_organ_onlyLable",date,".pdf"))
print(p2)
dev.off()


##################  umap_organ graph no legend
p2 <- ggplot(data=info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=0.1, size=0.00000001,aes(colour=organ))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_organ)+
	theme_classic()+
	labs(x='',y='',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank(),plot.title=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())+
	theme(legend.text=element_text(colour="black",size=8),legend.position="none",plot.title=element_text(hjust=0.5))
pdf(paste0("UMAP_organ_Nolegend_",date,".pdf"))
print(p2)
dev.off()
png(paste0("UMAP_organ_Nolegend_",date,".png"))
print(p2)
dev.off()



############################################ Fig2.1_Fig1B get umap_celltype graph
da <- data.frame(celltype=names(table(info$celltype)),class=c(paste0("0",1:9),10:length(table(info$celltype))),number=as.numeric(table(info$celltype)))
write.table(cbind(da,col_celltype),paste0("UMAP_celltype_color.csv"),col.names=T,row.names=F,sep=",",quote=F)
rownames(da) <- da[,1]
info$celltype_lable <- da[as.character(info$celltype),"class"]
#celltype lable
info %>%
  dplyr::group_by(celltype_lable) %>%
  summarize(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2)) -> centers
#info$celltype <- paste0(info$celltype,".",info$celltype_lable)

#### col_organ & celltype_lable
p1 <- ggplot(data=info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=0.1, size=0.00000001,aes(colour=organ))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_organ)+
	theme_classic()+
	labs(x='UMAP1',y='UMAP2',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank())+
	theme(legend.text=element_text(colour="black",size=12),legend.position="right",plot.title=element_text(hjust=0.5))+
  geom_point(data=centers, mapping=aes(x=UMAP_1, y=UMAP_2), size=0, alpha=0) +
  geom_text(data=centers, mapping=aes(label=celltype_lable), size=3, fontface="italic")

png(paste0("UMAP_organ_celltype.lable_",date,".png"),width=1000,height=500)
print(p1)
dev.off()

#### col_celltype & celltype_lable
p1 <- ggplot(data=info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=0.1, size=0.00000001,aes(colour=celltype))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_celltype)+
	theme_classic()+
	labs(x='UMAP1',y='UMAP2',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank())+
	theme(legend.text=element_text(colour="black",size=12),legend.position="right",plot.title=element_text(hjust=0.5))+
  geom_point(data=centers, mapping=aes(x=UMAP_1, y=UMAP_2), size=0, alpha=0) +
  geom_text(data=centers, mapping=aes(label=celltype_lable), size=3,fontface="italic")

pdf(paste0("UMAP_celltype_",date,".pdf"),width=30,height=8)
print(p1)
dev.off()
png(paste0("UMAP_celltype_",date,".png"),width=1800,height=500)
print(p1)
dev.off()


##################  umap_celltype graph only legend
only.info <- info[!duplicated(info$celltype),]
p1 <- ggplot(data=only.info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=1, size=6,aes(colour=celltype))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_celltype)+
	theme_classic()+
	labs(x='',y='',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),axis.line = element_blank())+
	theme(legend.text=element_text(colour="black",size=12),plot.title=element_text(hjust=0.5))+
  geom_point(data=centers, mapping=aes(x=UMAP_1, y=UMAP_2), size=0, alpha=0) +
  geom_text(data=centers, mapping=aes(label=celltype_lable), size=4, fontface="plain")

p2 <- p1+theme(legend.position="right")
pdf(paste0("UMAP_celltype_onlyLegend",date,".pdf"),width=25,height=8)
print(p2)
dev.off()
p2 <- p1+theme(legend.position="none")
pdf(paste0("UMAP_celltype_onlyLable",date,".pdf"))
print(p2)
dev.off()


##################  umap_celltype graph no legend
p2 <- ggplot(data=info,aes(x=UMAP_1, y=UMAP_2))+
	geom_point(alpha=0.1, size=0.00000001,aes(colour=celltype))+
	guides(colour=guide_legend(override.aes=list(size=3)))+
	scale_color_manual(values=col_celltype)+
	theme_classic()+
	labs(x='',y='',title="")+
	theme(plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches"),legend.title=element_blank(),plot.title=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks = element_blank(),axis.line = element_blank())+
	theme(legend.text=element_text(colour="black",size=8),legend.position="none",plot.title=element_text(hjust=0.5))
pdf(paste0("UMAP_celltype_Nolegend_",date,".pdf"))
print(p2)
dev.off()
png(paste0("UMAP_celltype_Nolegend_",date,".png"))
print(p2)
dev.off()



######################## Fig3.1_Fig1A get organ cellnum bar graph 
organ_cellnum <- as.data.frame(table(info$organ))
colnames(organ_cellnum) <- c("organ","Cell_Number")
#organ_cellnum <- organ_cellnum[rev(order(organ_cellnum$Cell_Number)),]
organ_cellnum$organ  <- factor(organ_cellnum$organ,as.character(organ_cellnum$organ))

pdf(paste0("organ_barplot_",date,".pdf"),width=22, height=10)
ggplot(organ_cellnum,aes(y=log10(Cell_Number),x=organ,fill=organ))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=col_organ)+
  theme_bw()+ labs(x='',y='Log10(Cell Number)',title="")+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"), axis.text.x=element_text(size=15,face="plain", vjust=0.5, hjust=1,angle=90), axis.text.y=element_text(size=15,face="plain"),axis.title=element_text(size=18,face="bold"))
dev.off()


######################## Fig3.2_Fig1C get celltype cellnum bar graph 
celltype_cellnum <- as.data.frame(table(info$celltype))
colnames(celltype_cellnum) <- c("celltype","Cell_Number")
#celltype_cellnum <- celltype_cellnum[rev(order(celltype_cellnum$Cell_Number)),]
celltype_cellnum$celltype  <- factor(celltype_cellnum$celltype,as.character(celltype_cellnum$celltype))

pdf(paste0("Celltype_barplot_",date,".pdf"),width=35, height=10)
ggplot(celltype_cellnum,aes(x=celltype,y=log10(Cell_Number),fill=celltype))+
  geom_bar(stat="identity",width=0.8)+
  scale_fill_manual(values=col_celltype)+
  theme_bw()+ labs(x='',y='Log10(Cell Number)',title="")+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"), axis.text.x=element_text(size=10,face="plain", vjust=0.5, hjust=1,angle=90), axis.text.y=element_text(size=10,angle=90,face="plain"),axis.title=element_text(size=12,face="bold"))
dev.off()

##################  celltype cellnum no legend
p <- ggplot(celltype_cellnum,aes(x=celltype,y=log10(Cell_Number),fill=celltype))+
  geom_bar(stat="identity",width=0.8)+
  theme_bw()+ labs(x='',y='Log10(Cell Number)',title="")+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"),axis.text.x=element_text(size=15,face="plain", vjust=0.5, hjust=1,angle=90), axis.text.y=element_text(size=15,angle=90,face="plain"),axis.title=element_text(size=25,face="bold") ,legend.position="none")
p1 <- p+scale_fill_manual(values=col_celltype)
pdf(paste0("Celltype_barplot_Nolegend_",date,".pdf"),width=20,height=8)
print(p1)
dev.off()
### grey
p2 <- p+ scale_fill_manual(values=rep("grey",length(col_celltype)))
pdf(paste0("Celltype_barplot_Nolegend_grey_",date,".pdf"),20,height=8)
print(p2)
dev.off()


######################## Fig3.3 get celltype bar graph of organ ratio
p <- ggplot(info,mapping=aes(x=celltype,fill=organ))+
  geom_bar(stat="count",position='fill',width=0.8)+
  scale_fill_manual(values=col_organ)+
  theme_bw()+
  labs(x='',y='Organ Ratio',title="")+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"),axis.text.x=element_text(colour="black",size=15,face="plain", vjust=0.5, hjust=1,angle=90), axis.text.y=element_text(colour="black",size=15,angle=90,face="plain"),axis.title=element_text(size=25,face="bold"))
pdf(paste0("Celltype_barplot_Organ.Ratio_",date,".pdf"),width=30, height=8)
print(p)
dev.off()
##################  organ ratio Nolegend
p1 <- p + theme(legend.position="none")
pdf(paste0("Celltype_barplot_Organ.Ratio_Nolegend_",date,".pdf"),width=20,height=8)
print(p1)
dev.off()



###################################################### Figure S1 ######################################################

################### Figure S1A: get boxplot of UMI and Gene count per organ
p <- ggplot(info,aes(x=factor(organ),y=log10(nUMI),fill=organ,color=organ))+
    geom_boxplot(alpha=1,width=0.7,outlier.size=0.1)+
    labs(x="",y="Log10(UMI Count)")+
    theme_classic()+
    scale_y_continuous(limits = c(2,5))+
    scale_fill_manual(values = col_organ)+
    scale_color_manual(values = rep("black",length(col_organ)))+scale_size_manual(values=0.05)+
    theme(axis.text.x=element_text(colour="black",size=15,angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(colour="black",size=15,angle=90),axis.title=element_text(colour="black",size=25,face="bold"))+ ylim(min(log10(info$nUMI))-0.005, max(log10(info$nUMI))+0.005) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))
pdf(paste0("Figure.S1_organ_boxplot_UMI_",date,".pdf"),width=23,height=8)
print(p)
dev.off()
pn <- p+theme(legend.position="none")
pdf(paste0("Figure.S1_organ_boxplot_UMI_Nolegend_",date,".pdf"),width=20,height=6.5)
print(pn)
dev.off()

p <- ggplot(info,aes(x=factor(organ),y=log10(nGene),fill=organ,color=organ))+
    geom_boxplot(alpha=1,width=0.7,outlier.size=0.1)+
    labs(x="",y="Log10(Gene Count)")+
    theme_classic()+
    scale_y_continuous(limits = c(2,5))+
    scale_fill_manual(values = col_organ)+
    scale_color_manual(values = rep("black",length(col_organ)))+scale_size_manual(values=0.05)+
    theme(axis.text.x=element_text(colour="black",size=15,angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(colour="black",size=15,angle=90),axis.title=element_text(colour="black",size=25,face="bold"))+ ylim(min(log10(info$nGene))-0.005, max(log10(info$nGene))+0.005) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))
pdf(paste0("Figure.S1_organ_boxplot_Gene_",date,".pdf"),width=23,height=8)
print(p)
dev.off()
pn <- p+theme(legend.position="none")
pdf(paste0("Figure.S1_organ_boxplot_Gene_Nolegend_",date,".pdf"),width=20,height=6.5)
print(pn)
dev.off()


################### Figure S1B: get boxplot of UMI and Gene count per celltype
p <- ggplot(info,aes(x=factor(celltype),y=log10(nUMI),fill=celltype,color=celltype))+
    geom_boxplot(alpha=1,width=0.8,outlier.size=0.1)+
    labs(x="",y="Log10(UMI Count)")+
    theme_classic()+
    scale_y_continuous(limits = c(2,5))+
    scale_fill_manual(values = col_celltype)+
    scale_color_manual(values = rep("black",length(col_celltype)))+scale_size_manual(values=0.05)+
    theme(axis.text.x=element_text(colour="black",size=15,angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(colour="black",size=15),axis.title=element_text(colour="black",size=25,face="bold"))+ ylim(min(log10(info$nUMI))-0.005, max(log10(info$nUMI))+0.005) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))
pdf(paste0("Figure.S1_celltype_boxplot_UMI_",date,".pdf"),width=30,height=8)
print(p)
dev.off()
pn <- p+theme(legend.position="none")
pdf(paste0("Figure.S1_celltype_boxplot_UMI_Nolegend",date,".pdf"),width=20,height=8)
print(pn)
dev.off()

p <- ggplot(info,aes(x=factor(celltype),y=log10(nGene),fill=celltype,color=celltype))+
    geom_boxplot(alpha=1,width=0.8,outlier.size=0.1)+
    labs(x="",y="Log10(Gene Count)")+
    theme_classic()+
    scale_y_continuous(limits = c(2,5))+
    scale_fill_manual(values = col_celltype)+
    scale_color_manual(values = rep("black",length(col_celltype)))+scale_size_manual(values=0.05)+
    theme(axis.text.x=element_text(colour="black",size=15,angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(colour="black",size=15),axis.title=element_text(colour="black",size=25,face="bold"))+ ylim(min(log10(info$nGene))-0.005, max(log10(info$nGene))+0.005) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))
pdf(paste0("Figure.S1_celltype_boxplot_Gene_",date,".pdf"),width=30,height=8)
print(p)
dev.off()
pn <- p+theme(legend.position="none")
pdf(paste0("Figure.S1_celltype_boxplot_Gene_Nolegend",date,".pdf"),width=20,height=8)
print(pn)
dev.off()



df <- info[,c("percent_mito","nUMI","nGene","celltype")]
df.mean <- aggregate(df[,1:3], by=list(df$celltype), FUN=mean)
df.median <- aggregate(df[,1:3], by=list(df$celltype), FUN=median)
out <- cbind(df.mean,df.median[,-1],as.numeric(table(list(df$celltype))))
colnames(out) <- c("celltype","mean_percent_mito","mean_nUMI","mean_nGene","median_percent_mito","median_nUMI","median_nGene","cell_number")
write.table(out,paste0("Table1.B_Glabal_celltype_",date,".xls",col.names=T,row.names=F,sep="\t",quote=F)

