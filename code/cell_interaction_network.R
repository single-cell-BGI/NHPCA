####################   ligand,receptor date

library(igraph)
library(RColorBrewer)

args <- commandArgs(T)
file <- args[1]
outpre <- args[2]
pros <- args[3]
date <- args[4]

pro1 <- unlist(strsplit(pros,split=","))[1]
pro2 <- unlist(strsplit(pros,split=","))[2]

#random_color <- c("#ffc87c","#E79076","#4682b4","#b699ff","#99B476","#c0d4f6","#6a534f","#3b3c36","#d53e4f","#9D6F90")
#random_color <- c(brewer.pal(12, "Paired"),brewer.pal(8, "Set2"))
random_color <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","9"="#D8A767","10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","19"="#E6C2DC","20"="#3D3D3D")

# read file
sig_mean_df <- read.table(file, sep='\t', header=T, as.is=T, check.names=F)
colnames(sig_mean_df) <- gsub(tail(unlist(strsplit(outpre,"/")),1)," ",colnames(sig_mean_df))
colnames(sig_mean_df)[grep("\\|",colnames(sig_mean_df))] <- gsub("_"," ",colnames(sig_mean_df)[grep("\\|",colnames(sig_mean_df))])
sig_mean_df <- sig_mean_df[setdiff(intersect(grep(pro1,sig_mean_df$interacting_pair),grep(pro2,sig_mean_df$interacting_pair)),grep("LGR4",sig_mean_df$interacting_pair)),]
#str(sig_mean_df)
rownames(sig_mean_df) <- sig_mean_df$interacting_pair
sig_mean_df <- sig_mean_df[13:ncol(sig_mean_df)]
#str(sig_mean_df)

sig_mean_df[is.na(sig_mean_df)] <- 0
sig_mean_df[sig_mean_df <= 0.13] <- 0
sig_mean_df[sig_mean_df > 0.13] <- 1
#str(sig_mean_df)
links <- c()
for(i in 1:nrow(sig_mean_df)){
	if(length(grep(pro1,sub("\\_.*","",rownames(sig_mean_df)[i])))>0){
	  link <- cbind(sub('\\|.*', '', colnames(sig_mean_df)[sig_mean_df[i,]>0]),sub('.*\\|', '', colnames(sig_mean_df)[sig_mean_df[i,]>0]))
	  links <- rbind(links,link)
	  }else{
	  link <- cbind(sub('.*\\|', '', colnames(sig_mean_df)[sig_mean_df[i,]>0]),sub('\\|.*', '', colnames(sig_mean_df)[sig_mean_df[i,]>0]))
	  links <- rbind(links,link)
	}
}
pair <- paste(links[,1],links[,2],sep="|")
if(length(grep("Low quality|Uncharacterized cell",pair))>0){
pair <- pair[-grep("Low quality|Uncharacterized cell",pair)]}else{
pair <- pair}
pair_num <- table(pair)

# set links
links <- data.frame(lig_rec=names(pair_num), pair.num=as.numeric(pair_num))
links$from <- sub('\\|.*', '', links$lig_rec)
links$to <- sub('.*\\|', '', links$lig_rec)
links$edge.color <- '#999999'
links <- links[c('from', 'to', 'pair.num')]
links$lig_rec <- NULL
links <- links[links$from != links$to,]
#str(links)

# set nodes
nodes <- data.frame(celltype=sort(unique(c(links$from,links$to))), stringsAsFactors=T)
rownames(nodes) <- nodes$celltype
nodes$vertex.color <- random_color[1:length(nodes$celltype)]
nodes$cell.num <- rep(100,nrow(nodes))

for(i in unique(links$from)) links[links$from==i, 'edge.color'] <- nodes[nodes$celltype==i, 'vertex.color']

# construct interaction net
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
pdf(paste0(outpre,pros,"_cell-cell_network_",date,".pdf"),width=20,height=20)
plot(net,
     vertex.color = V(net)$vertex.color,
     vertex.frame.color = "#444444",
     vertex.label.dist = V(net)$label.dist,
     vertex.label.degree = V(net)$lable.degree,
     vertex.label.cex = 2.5,
     vertex.label.color = "black",
     vertex.size = 40,
     edge.curved = 0.2,
     edge.color = E(net)$edge.color,
#     edge.width = (E(net)$pair.num)^3/(1.2),  ## 缩放线会变粗很多
     edge.width = (E(net)$pair.num)^2/(1.5),
     edge.arrow.size = 2,
     edge.arrow.size = ifelse((E(net)$pair.num)*(E(net)$pair.num)/5 >= rev(sort((E(net)$pair.num)*(E(net)$pair.num)/5))[3],6,2),
     edge.arrow.width = 1,
     layout=layout_in_circle)
dev.off()

### bigArrow
pdf(paste0(outpre,pros,"_cell-cell_network_bigArrow_",date,".pdf"),width=20,height=20)
plot(net,
     vertex.color = V(net)$vertex.color,
     vertex.frame.color = "#444444",
     vertex.label.dist = V(net)$label.dist,
     vertex.label.degree = V(net)$lable.degree,
     vertex.label.cex = 2.5,
     vertex.label.color = "black",
     vertex.size = 40,
     edge.curved = 0.2,
     edge.color = E(net)$edge.color,
     edge.width = (E(net)$pair.num)^2/(1.5),
     edge.arrow.size = 5,
     edge.arrow.size = ifelse((E(net)$pair.num)*(E(net)$pair.num)/5 >= rev(sort((E(net)$pair.num)*(E(net)$pair.num)/5))[3],6,2),
     edge.arrow.width = 1,
     layout=layout_in_circle)
dev.off()
