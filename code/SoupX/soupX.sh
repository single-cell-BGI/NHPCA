cat ../datasets_list | while read line
do
arr=(${line//\\t/ })
if [ -d "${arr[0]}" ]; then
        cd ${arr[0]}
else
        mkdir ${arr[0]} && cd ${arr[0]}
fi

name=${arr[0]}
output=${arr[1]}
raw="/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/4.bodyatlas/1.analysis/result/"${name}"/ALL"
flt=$(find /jdfssz1/ST_SUPERCELLS/P18Z10200N0350/Project/Monkey_Body_RNA_Response/0.HuaWei_5.0/Summary_Matrix | xargs ls -Rd| grep "${name}")

output="/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/1.Macaca_shanghai/huangbaoqian/07.soupX/08.no_remove/01.result_5.0"

filename=$name"_soupx.R"
cat>"${filename}"<<EOF
print(paste("Start time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
#.libPaths(c("/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/anaconda3/lib/R/library","/hwfssz1/ST_PRECISION/USER/zouxuanxuan/softwares/R/R-3.5.2/lib64/R/library","/hwfssz1/ST_PRECISION/PUB/softwares/R-3.5.1/lib64/R/library"))

.libPaths(c("/hwfssz1-tmp/ST_MCHRI/STEMCELL/USER/wuliang2/BackupData/wuliang2/biosoftware/anaconda3/lib/R/library"))

library(SoupX)
library(Seurat)

options(future.globals.maxSize = 100000 * 1024^3)
output <- paste("$output","$name",sep="/")
setwd(output)

toc <- Read10X("$flt",gene.column=1)
tod <- Read10X("$raw",gene.column=1)

tod <- tod[rownames(toc),]

all <- toc
all <- CreateSeuratObject(all)

all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)

matx <- all@meta.data

sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx\$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)
#samp <- c("Bone_marrow","PBMC","Pineal_gland")
if (length(grep("Bone_marrow|PBMC|Pineal_gland","$name"))==0) {
    if (unique(sc\$metaData\$rho) < 0.2) {
        sc = setContaminationFraction(sc, 0.2)
    }
}

out = adjustCounts(sc)
saveRDS(sc,"sc.rds")

DropletUtils:::write10xCounts("./soupX_matrix", out,version="3")

EOF

filename="r"$name"_soupx.sh"
cat>"${filename}"<<EOF
/hwfssz5/ST_PRECISION/OCG/wuliang2_SingleCell/RNA_pipeline/scRNA_V2.3/PISA-master/Rscript $name"_soupx.R"
#/hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/anaconda3/bin/Rscript $name"_soupx.R"
EOF

#qsub -cwd -binding linear:1 -l vf=1g,num_proc=1 -P P20Z10200N0059 -q st.q $filename
qsub -cwd -binding linear:1 -l vf=2g,num_proc=1 -P P20Z10200N0059 -q st_short.q $filename
cd $output

done
