###AnnotFile
for j in {1..22}
      do
        python make_annot.py --bed-file celltype.bed --bimfile 1000G.EUR.QC.${j}.bim --annot-file celltype.${j}.annot.gz
      done
###LDscore
for j in {1..22}
      do
       python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${j} --ld-wind-cm 1 --annot celltype.${j}.annot.gz --thin-annot --print-snps hapmap3_snps/hm.${j}.snp --out celltype.${j}
      done
