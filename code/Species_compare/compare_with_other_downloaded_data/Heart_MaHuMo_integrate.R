library(Seurat)
library(dplyr)
library(ggplot2)
organ <- "PBMC"
###########data for human#############



h_all<-readRDS("/hwfssz5/ST_PRECISION/TOMCAT/body_atlas/downloaddata/LGR5/result/heart_human.save.rds")
tmp<-h_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
h_all<-CreateSeuratObject(tmp)
h_all$Species <- rep(x="Human",times=nrow(h_all@meta.data))

h_all

rm(h1)
rm(h2)
rm(h1_o)
rm(h2_o)
###########data for macaca#############
ma_all<-readRDS("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/10.ANNO_RDS_0724/RDS/Anno_0724_Heart.rds")
ma_all$Species <- rep(x="Macaca",times=nrow(ma_all@meta.data))


###########data for mouse#############
mus_all<-readRDS("/hwfssz5/ST_PRECISION/TOMCAT/body_atlas/downloaddata/LGR5/result/heart_mouse.save.rds")
tmp<-mus_all@assays$RNA@counts
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
mus_all<-CreateSeuratObject(tmp)

mus_all$Species <- rep(x="Mouse",times=nrow(mus_all@meta.data))
rm(mus1)
mus_all
#############genes_selected#############

geneall<-c(rownames(h_all),rownames(ma_all),rownames(mus_all))
length(geneall)
a<-table(geneall)
keep<-names(a[a==3])
length(keep)
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


saveRDS(Brain.integrated.Find,file=paste0("Seurat_HuMaMo","heart",".integrated.Find.less.RDS"))











