library(JASPAR2018)

species <- 9606
opts <- list()
opts["species"] <- species
out  <- TFBSTools::getMatrixSet(JASPAR2018, opts)
names(out) <- paste(names(out), TFBSTools::name(out),sep = "_")
motif=out
seq1=read.table("seq1.txt",sep = "\t")
motif_ix_1 <- matchMotifs(motif, as.character(seq1[2,1])) 
r1=t(as.matrix(motifMatches(motif_ix_1)))
out=as.data.frame(cbind(r1)
out$TF=rownames(out)
write.table(out,"ACE2_peak_Motif.xls",sep = "\t",quote = FALSE,row.names = FALSE)
