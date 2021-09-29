args<-commandArgs(T)
library(Seurat)
library(dplyr)
library(ggplot2)
organ <- args[1]
###########data for human#############




all<-readRDS(args[2])
DefaultAssay(object = all) <- "RNA"
all$label<-paste0(all$Tissue,"_",all$Sample)
all@meta.data$donor2<-"batch1"
all@meta.data$donor2[all$label == "Visceral_adipose_MM2"]<-"batch2"

all@meta.data$donor2[all$label == "Diaphragm_MM2"]<-"batch2"
all@meta.data$donor2[all$label == "Bladder_FM3"]<-"batch2"


#all@meta.data$donor[which(all@meta.data$donor=="monkey2")]<-"monkey3"
obj.list<-SplitObject(all,split.by="donor2")
for (i in 1:length(x = obj.list)) {
    obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
reference.list <- obj.list
Brain.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Brain.integrated <- IntegrateData(anchorset = Brain.anchors, dims = 1:30)
library(ggplot2)
library(cowplot)

DefaultAssay(object = Brain.integrated) <- "integrated"
Brain.integrated <- ScaleData(object = Brain.integrated, verbose = FALSE)
Brain.integrated <- RunPCA(object = Brain.integrated, npcs = 30, verbose = FALSE)
Brain.integrated <- RunUMAP(object = Brain.integrated, reduction = "pca",dims = 1:30)
Brain.integrated = FindNeighbors(object = Brain.integrated,dims = 1:30)
Brain.integrated.Find <- FindClusters(object = Brain.integrated,resolution = 0.5)
DefaultAssay(Brain.integrated.Find)="RNA"
#Brain.integrated.Find@meta.data$mustype<-"fast muscle"
#Brain.integrated.Find@meta.data$mustype[which(Brain.integrated.Find@meta.data$seurat_clusters=="5")]<-"slow muscle"
#Brain.integrated.Find@meta.data$mustype[which(Brain.integrated.Find@meta.data$seurat_clusters=="6")]<-"slow muscle"
#Brain.integrated.Find@meta.data$mustype[which(Brain.integrated.Find@meta.data$seurat_clusters=="8")]<-"slow muscle"
#Brain.integrated.Find@meta.data$mustype[which(Brain.integrated.Find@meta.data$seurat_clusters=="9")]<-"slow muscle"

saveRDS(Brain.integrated.Find,file=paste0("Seurat_commoncell",organ,".integrated.Find.RDS"))
mycolor= c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#000011","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#616114","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8")


pdf(paste0("commoncell",organ,"_Cluster_Batch.less.pdf"),width=10,height=8)

DimPlot(object = Brain.integrated.Find, reduction = "umap",pt.size = 0.5,label=T,label.size=5)
DimPlot(object = Brain.integrated.Find, reduction = "umap",group.by = "organ",pt.size = 0.5,label=T,label.size=5)
DimPlot(object = Brain.integrated.Find, reduction = "umap",group.by = "donor2",pt.size = 0.5,label=F)
dev.off()

pdf(paste0("commoncell",organ,"_Cluster_BatchSplit.less.pdf"),width=21,height=8)
DimPlot(Brain.integrated.Find, reduction = "umap", split.by = "donor",pt.size = 0.5,no.legend = FALSE,label=T)
dev.off()

Idents(Brain.integrated.Find)<-Brain.integrated.Find$organ
DEG<-FindAllMarkers(Brain.integrated.Find,only.pos=T,logfc.threshold=0.1,min.pct=0.1)
write.table(DEG,file=paste0("commoncell_",organ,"_raw.DEG.xls"),sep="\t",quote=F)
DEG$FC<-exp(DEG$avg_logFC)
DEG<-DEG[DEG$FC>1.25,]
DEG<-DEG[DEG$pct.1>0.25,]
DEG<-DEG[DEG$p_val_adj<0.01,]
write.table(DEG,file=paste0("commoncell_",organ,"_cut.DEG.xls"),sep="\t",quote=F)

ave<-AverageExpression(Brain.integrated.Find,assays="RNA")
write.table(ave$RNA,file=paste0("commoncell_",organ,"_ave.xls"),sep="\t",quote=F)

Brain.integrated.Find<-Brain.integrated.Find.bak
Brain.integrated.Find<-subset(Brain.integrated.Find,idents=c("Diaphragm","Peritoneum","Tongue"))
Idents(Brain.integrated.Find)<-Brain.integrated.Find$mustype
DEG_slow<-FindAllMarkers(Brain.integrated.Find,min.pct = 0.25,logfc.threshold=0.5,only.pos = T)
write.table(DEG_slow,file=paste0("commoncell_",organ,"_slowfast.DEG.xls"),sep="\t",quote=F)

library(dplyr)
top30 <- DEG_slow %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DoHeatmap(rds, features = top30$gene)
dev.off()


Brain.integrated.Find$label<-paste0(Brain.integrated.Find$mustype,"_",Brain.integrated.Find$organ)
Brain.integrated.Find_slow<-subset(Brain.integrated.Find,subset=mustype=="slow muscle")
Brain.integrated.Find_fast<-subset(Brain.integrated.Find,subset=mustype=="fast muscle")
Idents(Brain.integrated.Find_slow)<-Brain.integrated.Find_slow$organ
DEG_slow_cut<-FindAllMarkers(Brain.integrated.Find_slow,min.pct = 0.25,logfc.threshold=0.5,only.pos = T)
Idents(Brain.integrated.Find_fast)<-Brain.integrated.Find_fast$organ
DEG_fast_cut<-FindAllMarkers(Brain.integrated.Find_fast,min.pct = 0.25,logfc.threshold=0.5,only.pos = T)


top10_s <- DEG_slow_cut %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_f <- DEG_fast_cut %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Brain.integrated.Find_slow, features = top10_s$gene)
DoHeatmap(Brain.integrated.Find_fast, features = top10_f$gene)


gene_s<-c("SPOCK3","DKK2","MAMDC2")
gene_f<-c("MEIS2","VAT1L","CPAMD8")
VlnPlot(Brain.integrated.Find_slow,features = gene_s,pt.size=0,ncol = 1)
VlnPlot(Brain.integrated.Find_fast,features = gene_f,pt.size=0,ncol = 1)

Idents(Brain.integrated.Find)<-Brain.integrated.Find$label

id_list<-unique(Brain.integrated.Find@meta.data$label)
library(cowplot)
avg.all.cells<- AverageExpression(Brain.integrated.Find,slot="data",assays = "RNA")
tmp<-as.data.frame(avg.all.cells$RNA)
list_1<-unique(colnames(tmp))
list_1<-sort(list_1)
tmp<-tmp[,list_1]
co<-cor(tmp,method="spearman")
cop<-cor(tmp,method="pearson")
library(pheatmap)
pheatmap(co,show_colnames = T,show_rownames = T,fontsize = 20,
             cluster_rows = T,cluster_cols = T,scale="none")
dev.off()
pheatmap(cop,show_colnames = T,show_rownames = T,fontsize = 20,
             cluster_rows = T,cluster_cols = T,scale="none")
dev.off()


DEG_label<-FindAllMarkers(Brain.integrated.Find,min.pct = 0.5,logfc.threshold=0.69315,only.pos = T)
pos<-which(rownames(tmp) %in% DEG_label$gene)
tmp_cut<-tmp[pos,]
cop_cut<-cor(tmp_cut,method="pearson")
co_cut<-cor(tmp_cut,method="spearman")
pheatmap(co_cut,show_colnames = T,show_rownames = T,fontsize = 20,
             cluster_rows = T,cluster_cols = T,scale="none")
pheatmap(cop_cut,show_colnames = T,show_rownames = T,fontsize = 20,
             cluster_rows = T,cluster_cols = T,scale="none")








