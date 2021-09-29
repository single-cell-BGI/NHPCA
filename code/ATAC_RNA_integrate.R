# Get the parameters
parser = argparse::ArgumentParser(description=cat("Script to integrate scATAC data and scRNA data\n\n# Authors:              Xi Dai\n# Contact information:  daixi@genomics.cn\n# Date:                 2020-07-27\n# R package version:    ArchR 0.9.5, Seurat 4.0.3\n\n"))
parser$add_argument('-A', '--ATACArchRProject', help='ATAC ArchR project')
parser$add_argument('-R', '--RNASeuratObject', help='RNA Seurat object, will use raw assay in Seurat object whose name is RNA')
parser$add_argument('-O', '--out', default='.', help='Out directory for results [default = .]')
parser$add_argument('-U', '--use', default='ArchR', help='ATAC and RNA integration methods to use, include ArchR, TransferData and IntegrateData [default = ArchR]')
args = parser$parse_args()


library(ArchR)
library(Seurat)
library(cowplot)


if(args$use == 'ArchR'){
	dir.create(paste0(args$out, "/ArchR"))
	seRNA <- readRDS(args$RNASeuratObject)
	DefaultAssay(seRNA) <- 'RNA'

	projHeme2 <- loadArchRProject(path = args$ATACArchRProject)
	
	projHeme2 <- addGeneIntegrationMatrix(
		ArchRProj = projHeme2, 
		useMatrix = "GeneScoreMatrix",
		matrixName = "GeneIntegrationMatrix",
		reducedDims = "IterativeLSI",
		seRNA = seRNA,
		addToArrow = FALSE,
		groupRNA = "Celltype",
		nameCell = "predictedCell_Un",
		nameGroup = "predictedGroup_Un",
		nameScore = "predictedScore_Un"
	)
	
	saveArchRProject(ArchRProj = projHeme2, outputDirectory = paste0(args$out, "/ArchR"), load = FALSE)


}else{
	pbmc_rna <- readRDS(args$RNASeuratObject)
	DefaultAssay(pbmc_rna) <- 'RNA'
	pbmc_rna <- FindVariableFeatures(object = pbmc_rna, nfeatures = 2000, verbose = FALSE)
	pbmc_rna <- ScaleData(object = pbmc_rna, verbose = FALSE)
	genesUse <- VariableFeatures(object = pbmc_rna)
	
	projHeme2 <- loadArchRProject(path = args$ATACArchRProject)
	
	GeneScoreMatrix <- getMatrixFromProject(
		ArchRProj = projHeme2,
		useMatrix = "GeneScoreMatrix"
	)
	gene_score <- assays(GeneScoreMatrix)$GeneScoreMatrix
	rownames(gene_score) <- rowData(GeneScoreMatrix)$name
	gene_score <- gene_score[rownames(gene_score) %in% genesUse,]
	
	mat <- log(gene_score + 1)
	pbmc <- CreateSeuratObject(
		counts = mat,
		assay = 'GeneScore',
		project = 'ATAC',
		min.cells = 1,
		meta.data = as.data.frame(projHeme2@cellColData)
	)
	pbmc <- ScaleData(pbmc, verbose = FALSE)
	
	UMAP <- getEmbedding(projHeme2, embedding = "UMAP", returnDF = TRUE)
	UMAP <- as.matrix(UMAP)
	pbmc[["umap"]] <- CreateDimReducObject(
		embeddings = UMAP, 
		key = "UMAP_", 
		assay = DefaultAssay(pbmc)
	)
	
	rD <- getReducedDims(
		ArchRProj = projHeme2, 
		reducedDims = "IterativeLSI",
		corCutOff = 1
	)
	pbmc[["lsi"]] <- CreateDimReducObject(
		embeddings = rD, 
		key = "LSI_", 
		assay = DefaultAssay(pbmc)
	)
	
	transfer.anchors <- FindTransferAnchors(
		reference = pbmc_rna,
		query = pbmc,
		reduction = 'cca'
	)
	predicted.labels <- TransferData(
		anchorset = transfer.anchors,
		refdata = pbmc_rna$Celltype,
		weight.reduction = pbmc[['lsi']]
	)
	pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
	
	p1 <- DimPlot(object = pbmc, label = TRUE, group.by = 'Clusters')
	p2 <- DimPlot(object = pbmc, label = TRUE, group.by = 'predicted.id')
	
	pdf(paste0(args$out, "/ATAC_umap.pdf"),width = 16,height = 8)
	plot_grid(p1, p2)
	dev.off()

	pdf(paste0(args$out, "/ATAC_umap.pdf"),width = 16,height = 8)
	plot_grid(p1, p2)
	dev.off()
	
	pbmc$Tech <- "ATAC"
	pbmc_rna$Tech <- "RNA"


	if(args$use == 'TransferData'){
		dir.create(paste0(args$out, "/TransferData"))
		refdata <- GetAssayData(pbmc_rna, assay = "RNA", slot = "data")[genesUse, ]
		imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
									weight.reduction = pbmc[["lsi"]])
		
		pbmc[["RNA"]] <- imputation
		DefaultAssay(pbmc) <- 'RNA'
		
		coembed <- merge(x = pbmc_rna, y = pbmc)
		coembed <- ScaleData(coembed, features = genesUse, do.scale = FALSE)
		coembed <- RunPCA(coembed, features = genesUse, verbose = FALSE)
		coembed <- RunUMAP(coembed, dims = 1:30)
		coembed <- RunTSNE(coembed, dims = 1:30)
		coembed$Merged_cluster <- ifelse(!is.na(coembed$Celltype), coembed$Celltype, coembed$predicted.id)
		
		p1 <- DimPlot(coembed, group.by = "Tech", reduction = "umap") 
		p2 <- DimPlot(coembed, group.by = "Merged_cluster", reduction = "umap", label = TRUE, repel = TRUE)
		p3 <- DimPlot(coembed, group.by = "Tech", reduction = "tsne") 
		p4 <- DimPlot(coembed, group.by = "Merged_cluster", reduction = "tsne", label = TRUE, repel = TRUE)
	
		pdf(paste0(args$out, "/TransferData/coembed_umap.pdf"), width = 16, height = 8)
		plot_grid(p1, p2)
		dev.off()
		
		pdf(paste0(args$out, "/TransferData/coembed_tsne.pdf"), width = 16, height = 8)
		plot_grid(p3, p4)
		dev.off()

		pdf(paste0(args$out, "/TransferData/coembed_umap.pdf"), width = 16, height = 8)
		plot_grid(p1, p2)
		dev.off()
		
		pdf(paste0(args$out, "/TransferData/coembed_tsne.pdf"), width = 16, height = 8)
		plot_grid(p3, p4)
		dev.off()
	
		save.image(paste0(args$out, "/TransferData/Transfer.Rdata"))
	}
	
	
	if(args$use == 'IntegrateData'){
		dir.create(paste0(args$out, "/IntegrateData"))
		DefaultAssay(pbmc) <- 'GeneScore'
		reference.list <- c(pbmc_rna, pbmc)
		names(reference.list) <- c("RNA", "ATAC")
		rna_atac.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
		rna_atac_integrated <- IntegrateData(anchorset = rna_atac.anchors, dims = 1:30)
		
		rna_atac_integrated <- ScaleData(object = rna_atac_integrated, verbose = F)
		rna_atac_integrated <- RunPCA(object = rna_atac_integrated, verbose = F)
		
		rna_atac_integrated <- FindNeighbors(object = rna_atac_integrated, dims = 1:30)
		rna_atac_integrated <- FindClusters(object = rna_atac_integrated, resolution = 0.5)
		
		rna_atac_integrated <- RunUMAP(object = rna_atac_integrated, reduction = "pca", dims = 1:30)
		rna_atac_integrated <- RunTSNE(object = rna_atac_integrated, reduction = "pca", dims = 1:30)
		
		rna_atac_integrated$Merged_cluster <- ifelse(!is.na(rna_atac_integrated$Celltype), rna_atac_integrated$Celltype, rna_atac_integrated$predicted.id)
		
		p1 <- DimPlot(rna_atac_integrated, group.by = "Tech", reduction = "umap")
		p2 <- DimPlot(rna_atac_integrated, group.by = "Merged_cluster", reduction = "umap", label = TRUE, repel = TRUE)
		p3 <- DimPlot(rna_atac_integrated, group.by = "Tech", reduction = "tsne")
		p4 <- DimPlot(rna_atac_integrated, group.by = "Merged_cluster", reduction = "tsne", label = TRUE, repel = TRUE)
		
		pdf(paste0(args$out, "/IntegrateData/integrated_umap.pdf"), width = 20, height = 10)
		plot_grid(p1, p2)
		dev.off()

		pdf(paste0(args$out, "/IntegrateData/integrated_umap.pdf"), width = 20, height = 10)
		plot_grid(p1, p2)
		dev.off()
	
		pdf(paste0(args$out, "/IntegrateData/integrated_tsne.pdf"), width = 20, height = 10)
		plot_grid(p3, p4)
		dev.off()

		pdf(paste0(args$out, "/IntegrateData/integrated_tsne.pdf"), width = 20, height = 10)
		plot_grid(p3, p4)
		dev.off()
	
		save.image(paste0(args$out, "/IntegrateData/Integrate.Rdata"))
	}
}
