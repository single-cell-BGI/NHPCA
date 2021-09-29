#!env python3
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import sys
import anndata.logging
import gc
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import random

random.seed(int(sys.argv[5]))

os.makedirs(sys.argv[1]) #create the output directory
os.chdir(sys.argv[1]) 
os.makedirs("write")

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
results_file = 'write/Organ45.h5ad'  # the file that will store the analysis results
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300, frameon=False)  # low dpi (dots per inch) yields small inline figures
# sc.settings.set_figure_params(dpi=80, facecolor='white')

mito_genes = ['ND1','ND2','ND3','ND4','ND4L','ND5','ND6','COX1','COX2','COX3','ATP6','ATP8','CYTB']

anndata.logging.print_memory_usage()

### input gene expression matrix
path = sys.argv[2] #原始矩阵目录
files = os.listdir(path) #得到文件夹下的所有文件名称
organ = [f[7:-4] for f in files] # extract organ name
adatas = []
for file in files: #遍历文件夹
     if not os.path.isdir(file): #判断是否是文件夹，不是文件夹才打开
          sub_adata = sc.read_csv(path+"/"+file, delimiter='\t', first_column_names=None, dtype='float32').T
          sub_adata.X = csr_matrix(sub_adata.X)
          print(f"{file[7:-4]}:{sub_adata.X.shape}")
          adatas.append(sub_adata) #每个文件的文本存到list中
adata = ad.concat(adatas, join="outer")
adata.obs['organ'] = ad.concat(adatas, label="organ", keys=organ).obs

del sub_adata
gc.collect()
anndata.logging.print_memory_usage()

### Preprocessing
#Show those genes that yield the highest fraction of counts in each single cells, across all cells.
sc.pl.highest_expr_genes(adata, n_top=20)
print(adata.X.shape)

sc.pp.filter_cells(adata, min_genes=1) # same for single organ 
sc.pp.filter_genes(adata, min_cells=3)
print(adata.X.shape)

### save raw count matrix of all cells
adata.write('write/Organ45_concatenate_rawcount_mat.h5ad',compression='gzip')

adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True, save='_QC_Organ45.pdf') ### bug
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_cm_Organ45.pdf')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_cg_Organ45.pdf')

adata

#adata = adata[adata.obs['n_genes'] < 4000, :]
#adata = adata[adata.obs['percent_mito'] < 0.3, :]

sc.pp.normalize_total(adata, target_sum=1e4)
### save normalize matrix of 43 organs
# adata.write('write/Organ45_concatenate_normalize_mat.h5ad',compression='gzip')
#pd.DataFrame(data=adata.X.T.todense(), index=adata.var_names,columns=adata.obs_names).to_csv('write/normalize_mat.tsv',sep="\t")
sc.pp.log1p(adata)
### save normalize and ln(x+1) matrix of 43 organs
adata.write('write/Organ45_concatenate_normalize_log_mat.h5ad',compression='gzip')
#pd.DataFrame(data=adata.X.T.todense(), index=adata.var_names,columns=adata.obs_names).to_csv('write/normalize_log_mat.tsv',sep="\t")
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)
### save scale matrix
adata.write('write/Organ45_concatenate_scale_mat.h5ad',compression='gzip')
#pd.DataFrame(data=adata.X.T, index=adata.var_names,columns=adata.obs_names).to_csv('scale_mat.csv',sep="\t")

### Principal component analysis
sc.tl.pca(adata, n_comps=50, random_state=int(sys.argv[3]), svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
adata

### Computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=int(sys.argv[3]))

### Embedding the neighborhood graph
sc.tl.umap(adata, random_state=int(sys.argv[5]))

### Clustering the neighborhood graph
sc.tl.louvain(adata,resolution=float(sys.argv[4]),random_state=int(sys.argv[3])) # set resolution
sc.tl.leiden(adata,resolution=float(sys.argv[4]),random_state=int(sys.argv[3]))
sc.pl.umap(adata, color=['louvain','leiden'],save='_leiden_'+sys.argv[5]+'.pdf')
sc.pl.umap(adata, color=['organ'],save='_Color_by_organ'+sys.argv[5]+'.pdf')
#adata.write(results_file, compression='gzip')

#####################################################  Marker Genes Feature Plot (UMAP)  #####################################################
########## 1. Subcutaneous adipose
marker = ["LPIN1","COL5A3","ACSL1","CXCL14","APOD","PECAM1","VWF","LYVE1","PTPRC","CD163","KIT","MS4A2","TAGLN","MYH11"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Subcutaneous.adipose.pdf')
########## 2. Visceral adipose
marker = ["LYVE1","TAGLN","MYH11","MS4A1","IL7R","CD34","MSLN","CFD"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Visceral.adipose.pdf')
########## 3. Adrenal gland
marker = ["MYLK","FLT1","PROX1","CYP17A1","CYP11B2","TOP2A","CHGA","TH","DCN","LAMA2","ITLN1","PKHD1L1","C1QA","APBB1IP"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Adrenal.gland.pdf')
########## 4. Artery
marker = ["RGS5","PDGFRB","FLT1","VWF","EMCN","PECAM1","CDH5","LYVE1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Artery.pdf')
########## 5. Bladder
marker = ["UPK3A","KRT5","TP63","KRT15","MKI67","TAGLN","MYH11","VIM","DCN","LAMA2","FLT1","PECAM1","VWF","LYVE1","PTPRC","KLRD1","CD3D","CD3G"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Bladder.pdf')
########## 6. Bone marrow
marker = ["GNLY","CCL5","GZMK","GZMB","GNLY","CCL5","GAMA","SH2D1A","LILRA4","JCHAIN","MS4A1","CD24","CD19","PF4","GP6","ITGA2B","TIMP3","NRGN","PPBP","NEGR1","TMEM154","TFRC","RHAG","KCNN4","GATA1","CD44","PRDX2","ITGA4","SLC40A1","PKIG","LEF1","MZB1","SLAMF7","VCAM1","LEPR","MGP","LPL","CXCL12","DCN","FBN1","S100A8","S100A9","FCER1A","LST1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Bone.marrow.pdf')
########## 7. Bronchus
marker = ["FOXJ1","CCDC78","TUBB1","TP73","AGER","MUC5B","SPDEF","SCGB3A2","PTPRC","CD163","SFTPB,SFTPC","LPO","KRT5","TP63","KRT15","FLT1","PECAM1","VWF","PROX1","DCN","LAMA2","ITLN1","MSLN","PKHD1L1","NEGR1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Bronchus.pdf')
########## 8. Carotid
marker = ["DCN","LAMA2","COL1A1","ATXN1","LPIN1","COL5A3","ACSL1","MYH11","TAGLN","FLT1","PECAM1","VWF","PROX1","CD3D","CD3G","IL7R","PDGFRB","KCNJ8","RGS5"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Carotid.pdf')
########## 9. Cerebellum
marker = ["PPP1R17","GABRA6","SLC17A7","SLC6A5","GRM2","SST","PRKCD","SORCS3","PTPRK","PRKCD","NXPH1","CDH22","EOMES","GDF10","AQP4","GJA1","SLC1A3","SLC1A2","MBP","MAG","MOBP","PLP1","PDGFRA","ITGAM","PTPRC","C1QB","CLDN5","PROM1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Cerebellum.pdf')
########## 10. Colon
marker = ["MZB1","SLAMF7","CD3D","CD3G","MS4A1","SDC1,KRT20,CDX2","MNDA","LYZ","PLP1","KIT","MS4A2","TOP2A"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Colon.pdf')
########## 11. Diaphragm
marker = ["NEB","TTN","LGR5","NEB","TTN","ITLN1","MSLN","PKHD1L1","DCN","LAMA2","COL1A1","FLT1","PECAM1","VWF","PROX1","LPIN1","COL5A3","ACSL1","PAX7","MYLK","DCN","LAMA2"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Diaphragm.pdf')
########## 12. Duodenum
marker = ["MZB1","SLAMF7","KRT19","EPCAM,","PLP1","PECAM1","VWF","LYVE1","MAP2","CD3D","CD3G","DCN","SLC1A2","SLC1A3","MNDA","LYZ"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Duodenum.pdf')
########## 13. Epididymis
marker = ["EPCAM","KRT19","KRT5","TP63","KRT15","MYH11","TAGLN","FLT1","PECAM1","VWF","PROX1","ATP6V1B1","CD3D","CD3G","IL7R","DCN","LAMA2","MYLK","LPIN1","COL5A3","ACSL1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Epididymis.pdf')
########## 14. Esophagus
marker = ["MYH3","MYH8","EPCAM","KRT19","KRT5","TP63","KRT15","CD3D","CD3G","IL7R","FLT1","PECAM1","VWF","PROX1","LPIN1","COL5A3","ACSL1","MYH11","TAGLN","MNDA","CD14","FCGR3A"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Esophagus.pdf')
########## 15. Fallopian tube
marker = ["PAX8","MKI67","MYLK","PKHD1L1","DCN","MYH11","TAGLN","LYVE1","PECAM1","GZMB","CCDC78"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Fallopian.tube.pdf')
########## 16. Gallbladder
marker = ["PECAM1","VWF","LYVE1","TAGLN","MYH11,","EPCAM","KRT19","PIGR","CYP2B6"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Gallbladder.pdf')
########## 17. Heart
marker = ["TNN2","MYH7","DCN","LAMA2","COL1A1","MYLK","DCN","LAMA2","FLT1","PECAM1","VWF","PROX1","LPIN1","COL5A3","ACSL1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Heart.pdf')
########## 18. Kidney
marker = ["PLA2R1","EGF","LRP2","FXYD4","FLT1","VWF","EMCN","SLC8A1","SLC12A3","FOXI1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Kidney.pdf')
########## 19. Liver
marker = ["SEC16B","ALB","KRT19","PDGFRB","DCN","FLT1","VWF","EMCN","PROM1","PTPRC","C1QA","CD163","CD14","MNDA","CD3D","CD3G","CCL5"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Liver.pdf')
########## 20. Lung
marker = ["FOXJ1","SCGB3A2","AGER","SFTPB","DCN","PTPRC","C1QA","FLT1","VWF","EMCN"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Lung.pdf')
########## 21. Lymphonode
marker = ["CD3D","CD3G","IL7R","DCN","LAMA2","MS4A1","CD24","CD19","FLT1","PECAM1","VWF","PROX1","MNDA","LYZ","LAMP3","LILRA4"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Lymphonode.pdf')
########## 22. Neocortex
marker = ["SLC17A7","GAD1","SLC1A3","SLC1A2","MBP","MAG","MOBP","PLP1","PDGFRA","ITGAM","PTPRC","C1QB","FLT1","CLDN5"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Neocortex.pdf')
########## 23. Nigra
marker = ["PDGFRA","TH","SLC6A3","ALDH1A1","SLC18A2","DDX","MBP","MAG","MOBP","PLP1","GAD1","GAD2"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Nigra.pdf')
########## 24. Ovary
marker = ["MKI67","NR5A2","TAGLN","CYP17A1","DCN","PKHD1L1","CCDC78","INHA"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Ovary.pdf')
########## 25. Pancreas
marker = ["GCG","INS","KRT19","DCN","FLT1","VWF","PECAM1","PROX1","PTPRC","C1QA","CD163"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Pancreas.pdf')
########## 26. PBMC
marker = ["MS4A1","GZMA","GZMK","GZMB","LEF1","PF4","FCGR3A"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_PBMC.pdf')
########## 27. Peritoneum
marker = ["NEB","TTN","DCN","LAMA2","COL1A1","LGR5","NEB","TTN","LPIN1","COL5A3","ACSL1","MYLK","DCN","LAMA2","FLT1","PECAM1","VWF","PROX1","PDGFRB","KCNJ8","RGS5","PAX7"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Peritoneum.pdf')
########## 28. Pigmentary_epithelium_choroid_plexus
marker = ["RPE65","BEST1","MLANA","PMEL","TYRP1","DCN","LAMA2","COL1A1","FLT1","PECAM1","VWF","PROX1","OPN1LW"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Pigmentary_epithelium_choroid_plexus.pdf')
########## 29. Pineal gland
marker = ["TPH1","AANAT","GFAP","PLP1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Pineal.gland.pdf')
########## 30. Pituitary gland
marker = ["FLT1","PECAM1","DCN","GHRHR","GH1","PRL","FSHB","GNRHR","POMC","TSHB","ANXA1","S100A6"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Pituitary.gland.pdf')
########## 31. Prostate gland
marker = ["EPCAM","KRT19","TP63","TAGLN","MYH11","DCN","LAMA2","SLC26A4","APT6V0D2","PTPRC","CD3D","CD3G","EPCAM","MKI67","FLT1","PECAM1","VWF","PROX1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Prostate.gland.pdf')
########## 32. Retina
marker = ["RBPMS","SLC17A6","OPN1LW","GNAT1","C1QL2","GAD1","VSX2","CALB1","APOE","RLBP1","CD163"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Retina.pdf')
########## 33. Salivary gland
marker = ["PRR27","LPO","MUC7","KLK1","ATP6V1B1","KRT7","DCN","LAMA2","AGR2","MTUS2","MYLK","FLT1","PECAM1","VWF","PROX1","CD3D","CD3G","SLC17A4","KRT5","KRT15","TP63","JCHAIN"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Salivary.gland.pdf')
########## 34. Spermaduct
marker = ["SEMG2","KRT5","KRT15","TP63","TAGLN","MYH11","SLC26A4","APT6V0D2","DCN","FLT1","COL5A3"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Spermaduct.pdf')
########## 35. Spleen
marker = ["GNLY","CCL5","GZMK","GZMB","GNLY","CCL5","LILRA4","JCHAIN","MS4A1","CD24","CD19","LEF1","MZB1","SLAMF7"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Spleen.pdf')
########## 36. Stomach
marker = ["PECAM1","VWF","LYVE1","KRT19","EPCAM,","PLP1","PDGFRA","SLC1A2","SLC1A3","DCN","TAGLN","MYH11","MZB1","SLAMF7","KIT","MS4A2","TOP2A","MNDA","LYZ","GNLY","CCL5","MS4A1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Stomach.pdf')
########## 37. Spinal cord
marker = ["MBP","MAG","MOBP","PLP1","CD163","FLT1","SYT1","NEFH","NEFM","MPZ","PMP22","PRX","SLC1A2","SLC1A3","PDGFRA"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Spinal.cord.pdf')
########## 38. Thyriod gland
marker = ["TG","TPO","IYD","FLT1","VWF","EMCN","PECAM1","DCN","ACTA1","COL5A3"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Thyriod.gland.pdf')
########## 39. Trachea
marker = ["FOXJ1","TUBB1","TP73","CCDC78","SCGB1A1","KRT5","KRT14","TAGLN","MYH11","DES","DCN","LPO","ILTN1","PKHD1L1","FLT1","VWF","PECAM1","LPIN1","COL5A3","ACSL1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Trachea.pdf')
########## 40. Tongue
marker = ["NEB","TTN","KRT19","EPCAM,","DCN","LAMA2","COL1A1","FLT1","PECAM1","VWF","PROX1","LPIN1","COL5A3","ACSL1","MYH11","TAGLN","PAX7","CD3D","CD3G","IL7R","MYLK","DCN","LAMA2"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Tongue.pdf')
########## 41. Myofibroblast
marker = ["MS4A1","CD3D","CD3G","KRT5","TP63","TAGLN","MYH11","MYH3","MYH8","MYLK","DCN","LILRA4","JCHAIN","MZB1","SLAMF7","FLT1","PECAM1","VWF","PROX1","DCN","LAMA2","SCGB3A1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Myofibroblast.pdf')
########## 42. Uterus
marker = ["FLT1","PECAM1","EPCAM","KRT5","KRT14","KRT17","ITLN1","MSLN","PKHD1L1","TAGLN","DCN","LAMA2","GZMB","MS4A1"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Uterus.pdf')
########## 43. Vagina
marker = ["TAGLN","MYH11","NEB","FLT1","PECAM1","KRT13","KRT4"]
marker1 = [gene for gene in marker if (gene in adata.raw.var._stat_axis.values.tolist())]
sc.pl.umap(adata, color=marker1, sort_order=True, size=2, save='_Vagina.pdf')


### Get cluster, umap coordinates and other information of all cells
cluster=pd.concat([adata.obs['organ'][:],adata.obs['leiden'][:],adata.obs['louvain'][:],adata.obs['percent_mito'][:],adata.obs['n_counts'][:],adata.obs['n_genes'][:]],axis=1)
cluster["UMAP_1"]=adata.obsm["X_umap"][:,0]
cluster["UMAP_2"]=adata.obsm["X_umap"][:,1]
cluster.columns = ["organ","leiden","louvain","percent_mito","nUMI","nGene","UMAP_1","UMAP_2"]
cluster = cluster[["organ","leiden","louvain","UMAP_1","UMAP_2","percent_mito","nUMI","nGene"]]
cluster.to_csv('write/Scanpy_cluster_info.csv')

### save cluster color
#pd.DataFrame(adata.uns['louvain_colors']).to_csv("write/Scanpy_louvain_cluster_color.csv",header=None)
#pd.DataFrame(adata.uns['leiden_colors']).to_csv("write/Scanpy_leiden_cluster_color.csv",header=None)

### Finding marker genes
sc.settings.verbosity = 2  # reduce the verbosity

#sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save='.pdf')

### Get a table with the scores and groups.
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
## top25_ranked_genes_per_cluster 
#trgpc = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores', 'logfoldchanges','pvals_adj']})
#print(trgpc)
#trgpc.to_csv(path_or_buf="write/ranked_genes_per_leiden_cluster.csv")
trgpc = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores', 'logfoldchanges','pvals_adj']}).head(25)
print(trgpc)
trgpc.to_csv(path_or_buf="write/top25_ranked_genes_per_leiden_cluster.csv")

adata
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

anndata.logging.print_memory_usage()
