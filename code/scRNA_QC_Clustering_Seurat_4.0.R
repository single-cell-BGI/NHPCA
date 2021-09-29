### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input directory')
parser$add_argument('-D','--id', help='sample ID')
parser$add_argument('-G','--filterGene', help='gene number threhold')
parser$add_argument('-M','--filterMito', help='mito ratio threhold')
parser$add_argument('-O','--out', help='out directory')
args = parser$parse_args()

###Combine mtx
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)

setwd(args$input)

file_list=list.files(pattern = args$id)
file_list=paste0(file_list,"/soupX_matrix")

obj.list=list()

for(i in file_list){

input=Read10X(data.dir = i,gene.column = 1)
tag=i
obj=CreateSeuratObject(counts = input, min.cells = 3, min.features = 200)
obj$Batch=tag
mito_genes = c('ND5','COX3','ATP8','COX1','ND6','ND3','ND4L','COX2','ND1','CYTB','ATP6','ND4','ND2')
mito_genes_1 = intersect(mito_genes,rownames(obj))
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito_genes_1)
obj.list[[i]]=subset(obj, subset = nFeature_RNA > as.numeric(args$filterGene) & percent.mt < as.numeric(args$filterMito))
#obj.list[[i]]=subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
                    }

###QC & Clustering

#########
Find_doublet <- function(data){
        sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        nExp_poi <- round(0.05*ncol(data))
        p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
        data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
        colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
        data
}
#########

doublet=list()

for (i in 1:length(x = obj.list)) {
    obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
     obj.list[[i]] <- ScaleData(obj.list[[i]])
     obj.list[[i]] <- RunPCA(obj.list[[i]])
     obj.list[[i]] <- RunUMAP(obj.list[[i]], dims = 1:10)
     obj.list[[i]] <- Find_doublet(obj.list[[i]])
     doublet[[i]] <- obj.list[[i]]@meta.data[,-6]
     obj.list[[i]] <- subset(obj.list[[i]],subset=doublet_info=="Singlet")
}

doublet_df=do.call(rbind,doublet)
write.table(doublet_df,file = paste0(args$out,"/","Doublet_info_",args$id,".txt"),sep="\t",quote=FALSE)

reference.list <- obj.list
features <- SelectIntegrationFeatures(object.list = reference.list)
obj.anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = features)
obj.integrated <- IntegrateData(anchorset = obj.anchors)

library(ggplot2)
library(cowplot)

DefaultAssay(object = obj.integrated) <- "integrated"
obj.integrated <- ScaleData(object = obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(object = obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca",
    dims = 1:30)
obj.integrated = FindNeighbors(object = obj.integrated,k.param=40,dims = 1:30)
obj.integrated.Find <- FindClusters(object = obj.integrated,resolution = 0.5)
obj.integrated.Find <- RunTSNE(object = obj.integrated.Find,dims = 1:30,seed.use = 2)
DefaultAssay(obj.integrated.Find)="RNA"

saveRDS(obj.integrated.Find,file=paste0(args$out,"/","Seurat_",args$id,".integrated.Find.rds"))

inputmarkers <- FindAllMarkers(obj.integrated.Find)
write.table(inputmarkers,paste0(args$out,"/","Marker_",args$id,".xls"),sep="\t",quote = FALSE)

#count_F=as.data.frame(obj.integrated.Find@assays$RNA@counts)
#write.table(count_F,paste0(args$out,"/","Filtered_count_",args$id,".txt"),sep="\t",quote = FALSE)
