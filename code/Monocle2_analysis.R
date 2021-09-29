library(Seurat)
library(monocle)
library(viridis)

rds=readRDS("Anno_0724_Tissue.rds")
Idents(rds) <- rds$Celltype
rds <- subset(rds, idents = c("celltype1","celltype2"))
data <- as(as.matrix(rds@assays$RNA@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = rds@meta.data)
fd0 <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fd0)
cds <- newCellDataSet(data,
                      phenoData  = pd,
                      featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~Celltype",
                                      reducedModelFormulaStr = "~Batch",
                                      cores=8)
ordering.genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:500]

cds_of <- setOrderingFilter(cds, ordering.genes)
cds_of <- reduceDimension(cds_of, max_components = 2, method = 'DDRTree')
cds_of <- orderCells(cds_of)
saveRDS(cds_of,file="cds_of_Tissue_2.RDS")
