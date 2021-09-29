library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
organ <- "Brain_V2"
###########data for human#############

#setwd("/hwfssz1-tmp/ST_PRECISION/USER/zhuangzhenkun/F16ZQSB1SY2825/USER_bak/zhuangzhenkun/body_atlas/10.integrate/evolution/Three_species/data/brain_human_mtx")
#file_list=list.files(pattern = "mtx")
#data=list()
#for(i in file_list){
#        input=as.data.frame(fread(i))
#        tag=gsub("_mat_ortholog_Ma_gene.txt","",i)
#        colnames(input)=c("ID",paste(tag,colnames(input)[2:ncol(input)],sep = ":"))
#        data[[i]]=input
#}
#
#df1 <- Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ID"),data)
#df1$ID <- gsub("_","-",df1$ID)
#df1[is.na(df1)]=0
#
#dim(df1)
#rownames(df1)<-df1[,1]
#df1<-df1[,-1]
a<-readRDS("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/0.hwfssz1/Monkey_Atlas_RNA/3.Species/ANNO_Seurat_HuMaMoBrain_new.integrated.Find.RDS")
a<-subset(a,subset=Species == "Human")
a
df1<-a@assays$RNA@counts
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/1.script")
h_all<-CreateSeuratObject(df1)
h_all[["percent.mt"]] <- PercentageFeatureSet(h_all, pattern = "^MT-")
h_all <- subset(h_all, subset = nFeature_RNA > 500 & percent.mt < 10)
tmp<-h_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
#a<-sample(x=c(1:ncol(tmp)),size = 10000)
#tmp<-tmp[,a]
h_all<-CreateSeuratObject(tmp)
h_all$Species <- rep(x="Human",times=nrow(h_all@meta.data))

h_all

rm(df1)
###########data for macaca#############
macaca1<-readRDS("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/10.ANNO_RDS_0724/RDS/Anno_0724_Neocortex.rds")
ma_all<-CreateSeuratObject(macaca1@assays$RNA@counts)
mito_genes = c('ND5','COX3','ATP8','COX1','ND6','ND3','ND4L','COX2','ND1','CYTB','ATP6','ND4','ND2')
ma_all[["percent.mt"]] <- PercentageFeatureSet(ma_all, features = mito_genes)
ma_all <- subset(ma_all, subset = nFeature_RNA > 500 & percent.mt < 10)
tmp<-ma_all@assays$RNA@counts
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Crab.eating.macaque.gene.name,]
a<-sample(x=c(1:ncol(tmp)),size = 10000)
tmp<-tmp[,a]
ma_all<-CreateSeuratObject(tmp)
ma_all$Species <- rep(x="Macaca",times=nrow(ma_all@meta.data))
rm(macaca1)


ma_all
###########data for mouse#############
mus1<-read.table("/hwfssz5/ST_PRECISION/TOMCAT/body_atlas/downloaddata/brain/mouse_BGI/Mouse_Brain_BGI_mtx.tsv",header=T,row.names=1,sep="\t")
mus_all<-CreateSeuratObject(mus1)
mus_all[["percent.mt"]] <- PercentageFeatureSet(mus_all, pattern = "^Mt-")
mus_all <- subset(mus_all, subset = nFeature_RNA > 500 & percent.mt < 10)

tmp<-mus_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
#a<-sample(x=c(1:ncol(tmp)),size = 10000)
#tmp<-tmp[,a]
mus_all<-CreateSeuratObject(tmp)

mus_all$Species <- rep(x="Mouse",times=nrow(mus_all@meta.data))
rm(mus1)
mus_all
#############genes_selected#############

geneall<-c(rownames(h_all),rownames(ma_all),rownames(mus_all))
length(geneall)
a<-table(geneall)
keep<-names(a[a==3])
all<-merge(x=h_all,y=c(ma_all,mus_all))
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
Brain.integrated.Find <- FindClusters(object = Brain.integrated,resolution = 0.5)
DefaultAssay(Brain.integrated.Find)="RNA"


saveRDS(Brain.integrated.Find,file=paste0("Seurat_HuMaMo",organ,".new.integrated.Find.RDS"))
mycolor= c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#000011","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#616114","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8")


pdf(paste0("HuMaMo",organ,"_new_Cluster_Batch.pdf"),width=10,height=8)

DimPlot(object = Brain.integrated.Find, reduction = "umap",pt.size = 0.5,cols=mycolor, no.legend = FALSE,label=T,lable.size=5)
DimPlot(object = Brain.integrated.Find, reduction = "umap",group.by = "Species",pt.size = 0.5,cols=mycolor,no.legend = FALSE,label=F)
dev.off()

pdf(paste0("HuMaMo",organ,"_new_Cluster_BatchSplit.pdf"),width=21,height=8)
DimPlot(Brain.integrated.Find, reduction = "umap", split.by = "Species",pt.size = 0.5,no.legend = FALSE,label=T)
dev.off()











