library(Seurat)
library(dplyr)
library(ggplot2)
organ <- "Kidney_7"
###########data for human#############



h1<-readRDS("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/zhuangzhenkun/USER_bak/zhuangzhenkun/body_atlas/4.human_kidney_clustering/Human200.RDS")
h_all<-CreateSeuratObject(h1@assays$RNA@counts)
h_all
h_all[["percent.mt"]] <- PercentageFeatureSet(h_all, pattern = "^MT-")
h_all <- subset(h_all, subset = nFeature_RNA > 200 & percent.mt < 10)
tmp<-h_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
h_all<-CreateSeuratObject(tmp)
h_all$Species <- rep(x="Human",times=nrow(h_all@meta.data))
h_all

h2<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/4.guoguoji_HCL_MCA/HCL/RDS/Kidney.origin.seurat.RDS")
h_all2<-CreateSeuratObject(h2@assays$RNA@counts)
h_all2[["percent.mt"]] <- PercentageFeatureSet(h_all2, pattern = "^MT-")
h_all2 <- subset(h_all2, subset = nFeature_RNA > 200 & percent.mt < 10)
tmp<-h_all2@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
a<-sample(x=c(1:ncol(tmp)),size = 10000)
tmp<-tmp[,a]
h_all2<-CreateSeuratObject(tmp)
h_all2$Species <- rep(x="Human2",times=nrow(h_all2@meta.data))
h_all2

h3<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/7.final_figure/4.figure4.new/script/KidneyHuman.rds")
h_all3<-CreateSeuratObject(h3@assays$RNA@counts)
h_all3[["percent.mt"]] <- PercentageFeatureSet(h_all3, pattern = "^MT-")
h_all3 <- subset(h_all3, subset = nFeature_RNA > 200 & percent.mt < 10)
tmp<-h_all3@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
a<-sample(x=c(1:ncol(tmp)),size = 10000)
tmp<-tmp[,a]
h_all3<-CreateSeuratObject(tmp)
h_all3$Species <- rep(x="Human3",times=nrow(h_all3@meta.data))
h_all3



###########data for macaca#############
macaca1<-readRDS("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/10.ANNO_RDS_0724/RDS/Anno_0724_Kidney.rds")
ma_all<-CreateSeuratObject(macaca1@assays$RNA@counts)
mito_genes = c('ND5','COX3','ATP8','COX1','ND6','ND3','ND4L','COX2','ND1','CYTB','ATP6','ND4','ND2')
ma_all[["percent.mt"]] <- PercentageFeatureSet(ma_all, features = mito_genes)
ma_all <- subset(ma_all, subset = nFeature_RNA > 200 & percent.mt < 10)
tmp<-ma_all@assays$RNA@counts
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Crab.eating.macaque.gene.name,]
a<-sample(x=c(1:ncol(tmp)),size = 10000)
tmp<-tmp[,a]
ma_all<-CreateSeuratObject(tmp)
ma_all$Species <- rep(x="Macaca",times=nrow(ma_all@meta.data))
rm(macaca1)


ma_all
###########data for mouse#############
mus_all<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/7.final_figure/4.figure4.new/script/KidneyMouse.rds")
#mus_all<-CreateSeuratObject(mus1)
mus_all[["percent.mt"]] <- PercentageFeatureSet(mus_all, pattern = "^Mt-")
mus_all <- subset(mus_all, subset = nFeature_RNA > 200 & percent.mt < 10)

tmp<-mus_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
a<-sample(x=c(1:ncol(tmp)),size = 10000)
tmp<-tmp[,a]
mus_all<-CreateSeuratObject(tmp)

mus_all$Species <- rep(x="Mouse",times=nrow(mus_all@meta.data))
rm(mus1)
mus_all


m2<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/4.guoguoji_HCL_MCA/MCA/RDS/Kidney.origin.seurat.RDS")
mus_all2<-CreateSeuratObject(m2@assays$RNA@counts)
#mus_all<-CreateSeuratObject(mus1)
mus_all2[["percent.mt"]] <- PercentageFeatureSet(mus_all2, pattern = "^Mt-")
mus_all2 <- subset(mus_all2, subset = nFeature_RNA > 200 & percent.mt < 10)

tmp<-mus_all2@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
mus_all2<-CreateSeuratObject(tmp)
mus_all2$Species <- rep(x="Mouse2",times=nrow(mus_all2@meta.data))
mus_all2


m3<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/5.Tabula/TB.kidney.10X.RDS")
mus_all3<-CreateSeuratObject(m3)
#mus_all<-CreateSeuratObject(mus1)
mus_all3[["percent.mt"]] <- PercentageFeatureSet(mus_all3, pattern = "^Mt-")
mus_all3 <- subset(mus_all3, subset = nFeature_RNA > 200 & percent.mt < 10)
tmp<-mus_all3@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1/ST_SUPERCELLS/P20Z10200N0059/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
mus_all3<-CreateSeuratObject(tmp)
mus_all3$Species <- rep(x="Mouse3",times=nrow(mus_all3@meta.data))
mus_all3

#############genes_selected#############

geneall<-c(rownames(h_all),rownames(h_all2),rownames(h_all3),rownames(ma_all),rownames(mus_all),rownames(mus_all2),rownames(mus_all3))
length(geneall)
a<-table(geneall)
keep<-names(a[a==7])
length(keep)
all<-merge(x=h_all,y=c(h_all2,h_all3,ma_all,mus_all,mus_all2,mus_all3))
all
all<-subset(all,features=keep)
all
obj.list<-SplitObject(all,split.by="Species")
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
Brain.integrated = FindNeighbors(object = Brain.integrated,k.param=40,dims = 1:30)
Brain.integrated <- FindClusters(object = Brain.integrated,resolution = 0.5)
DefaultAssay(Brain.integrated)="RNA"


saveRDS(Brain.integrated,file=paste0("Seurat_HuMaMo",organ,".integrated.Find.RDS"))
mycolor= c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#000011","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#616114","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8")
DEG<-FindAllMarkers(Brain.integrated,only.pos=T)
write.table(DEG,file=paste0("Seurat_HuMaMo",organ,".integrated.Find.DEG.xls"),quote=F,sep="\t")

pdf(paste0("HuMaMo",organ,"_Cluster_Batch.pdf"),width=10,height=8)

DimPlot(object = Brain.integrated, reduction = "umap",pt.size = 0.5,cols=mycolor,label=T,label.size=5)
DimPlot(object = Brain.integrated, reduction = "umap",group.by = "Species",pt.size = 0.5,cols=mycolor,label=F)
dev.off()

pdf(paste0("HuMaMo",organ,"_Cluster_BatchSplit.pdf"),width=21,height=8)
DimPlot(Brain.integrated.Find, reduction = "umap", split.by = "Species",pt.size = 0.5,label=T)
dev.off()






kidney$Celltype<-"none"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "0"]<-"Proximal tubule"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "1"]<-"Ascending LOH"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "2"]<-"Descending LOH"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "3"]<-"Principal"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "4"]<-"Proximal tubule"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "5"]<-"DCTC"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "6"]<-"Intercalated"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "7"]<-"Endothelial"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "8"]<-"Descending LOH"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "9"]<-"Macrophage"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "10"]<-"Podocyte"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "11"]<-"Stromal"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "12"]<-"NKT"
kidney$Celltype[kidney$integrated_snn_res.0.5 == "13"]<-"Principal like"
kidney$Celltype[kidney$integrated_snn_res.1 == "9"]<-"CTC"

#Idents(kidney)<-kidney$Celltype
#kidney_cut<-subset(kidney,subset=Celltype == args[1])
#DE_list<-list()
#Idents(kidney)<-kidney$Species
#sp<-c("Human","Human2","Mouse","Mouse2")
#for (i in c(1:4)){
#	DEG<-FindMarkers(kidney_cut,ident.1=sp[i],ident.2="Macaca",min.pct=0.1,logfc.threshold=0.1)
#	DE_list[[i]]<-DEG
#	write.table(DEG,file=paste0(sp[i],"_vs_Macaca.DEG.xls"),quote=F,sep="\t")
#}
#saveRDS(DE_list,file=paste0(args[1],"_DE.list.RDS"))






