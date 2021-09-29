args<-commandArgs(T)
library(Seurat)
library(dplyr)
library(ggplot2)
organ <- args[1]
organ2<-args[2]
out<-args[3]
###########data for human#############
a<-readRDS(paste0("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/4.guoguoji_HCL_MCA/HCL/RDS/",organ,".origin.seurat.RDS"))
tmp<-a@assays$RNA@counts
rm(a)
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/HuMa_one2one_16758genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
h_all<-CreateSeuratObject(tmp)
h_all$Species <- rep(x="Human",times=nrow(h_all@meta.data))
h_all$Abbreviation <- "zHuman"
h_all
rm(h1)
rm(h2)
rm(h1_o)
rm(h2_o)
###########data for macaca#############
ma_all<-readRDS(paste0("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/10.ANNO_RDS_0724/RDS/Anno_0724_",organ2,".rds"))
ma_all$Species <- rep(x="Macaca",times=nrow(ma_all@meta.data))


###########data for mouse#############
a<-readRDS(paste0("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/6.three_species_integrate/4.guoguoji_HCL_MCA/MCA/RDS/",organ,".origin.seurat.RDS"))
tmp<-a@assays$RNA@counts
rm(a)
HuMa_ortholog<-read.table("/hwfssz1-tmp/ST_PRECISION/USER/pantaotao/projects/00.Macaca_Fascicularis_Atlas/0.refdata/0.orthologous_genes/MoMa_one2one_14717genename_list.txt",header=T,sep="\t")
HuMa_ortholog<-as.data.frame(HuMa_ortholog)
tmp<-tmp[rownames(tmp) %in% HuMa_ortholog$Gene.name,]
rownames(HuMa_ortholog)<-HuMa_ortholog[,1]
HuMa_ortholog1<-HuMa_ortholog[rownames(tmp),]
HuMa_ortholog1<-as.data.frame(HuMa_ortholog1)
rownames(tmp)<-HuMa_ortholog1$Crab.eating.macaque.gene.name
mus_all<-CreateSeuratObject(tmp)

mus_all$Species <- rep(x="Mouse",times=nrow(mus_all@meta.data))
mus_all$Abbreviation <- "zMouse"
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
Brain.integrated <- FindClusters(object = Brain.integrated,resolution = 0.5)
DefaultAssay(Brain.integrated)="RNA"
saveRDS(Brain.integrated,file=paste0(out,"/Seurat_HuMaMo",organ2,".integrated.Find.singlecell.guoguiji.RDS"))
pdf(paste0(out,"/Seurat_HuMaMo",organ2,".integrated.anno.pdf"),width=12,height=10)
DimPlot(Brain.integrated,group.by = "Abbreviation",pt.size=1)
DimPlot(Brain.integrated,pt.size=1)
DimPlot(Brain.integrated,split.by = "Species")
dev.off()


DEG<-FindAllMarkers(Brain.integrated,min.pct =0.2,only.pos=T,logfc.threshold = 0.5,return.thresh =0.05)
write.table(DEG,file=paste0(out,"/Seurat_HuMaMo",organ2,"singlecell.integrated.deg.xls"),)




#saveRDS(Brain.integrated.Find,file=paste0("Seurat_HuMaMo",organ2,".integrated.Find.guoguiji.RDS"))











