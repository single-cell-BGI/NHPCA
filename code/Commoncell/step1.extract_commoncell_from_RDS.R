args<-commandArgs(T)
library(Seurat)
setwd("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/10.ANNO_RDS_0724/RDS")
#rds_list=list.files(pattern = "integrated.Find.RDS")
#rds_list=c("Seurat_anno_Subcutaneous_adipose.integrated.Find.RDS","Seurat_anno_Visceral_adipose.integrated.Find.RDS")
ob.list<-list()
ref<-read.table(args[1],header=T,sep="\t")

organ<-unique(as.vector(ref$Tissue))

find_cell <- function(organ){
  ref_cut<-ref[ref$Tissue == organ,]
  ref_cut$Abbreviation<-as.vector(ref_cut$Abbreviation)
  rds_p<-paste0("Anno_0724_",organ,".rds")
  rds<-readRDS(rds_p)
  rds_p<-gsub(pattern = "Anno_0721_|.rds",replacement = "",x=rds_p)
  rds@meta.data$organ<-rep(x=rds_p,times=nrow(rds@meta.data))
  tmp<-rds@meta.data
  tmp<-tmp[as.vector(tmp$Abbreviation) %in% ref_cut$Abbreviation,]
  if(nrow(tmp) >= 200){
   tag=1
  rds_cut<-subset(rds,cell=rownames(tmp))
  }
  else{
   tag=0
   rds_cut<-rds
  }
  out<-list(object=rds_cut,ta=tag)
  return(out)
}
j=1
for(i in c(1:length(organ))){
  aa<-find_cell(organ[i])
  if(aa[[2]] == 1){
   ob.list[[j]] = aa[[1]]
   j=j+1
   }
}

ob.list1<-ob.list[[1]]
ob.list[[1]]<-NULL
ob.listo<-ob.list
merge_rds<-merge(x=ob.list1,y=ob.listo)
setwd("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/7.final_figure/2.figure2.new/1.commoncell")
write.table(merge_rds@meta.data,file=paste0(args[2],".all.metadata.xls"),quote=F,sep="\t")
saveRDS(merge_rds,file=paste0(args[2],".all.rds"))


#table(merge_rds$Batch)
#table(merge_rds$organ)
Brain.integrated<-merge_rds
rm(merge_rds)
Brain.integrated <- NormalizeData(object = Brain.integrated, verbose = FALSE)
Brain.integrated <- FindVariableFeatures(object = Brain.integrated)
Brain.integrated <- ScaleData(object = Brain.integrated, verbose = FALSE)
Brain.integrated <- RunPCA(object = Brain.integrated, npcs = 30, verbose = FALSE)
Brain.integrated <- RunUMAP(object = Brain.integrated, reduction = "pca",dims = 1:30)
Brain.integrated = FindNeighbors(object = Brain.integrated,dims = 1:30)
Brain.integrated <- FindClusters(object = Brain.integrated,resolution = 0.5)
saveRDS(Brain.integrated,file=paste0("Seurat_commoncell",args[2],".no_integrated.Find.RDS"))





