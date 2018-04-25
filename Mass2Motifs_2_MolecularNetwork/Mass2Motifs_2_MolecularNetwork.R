###########################################################
#                                                         #
# Match Mass2Motifs to mass spectral molecular networks   #
#                                                         #
###########################################################
##############################################
#  map Mass2Motifs on edges (multiple edges) #
##############################################

library(plyr)

edges <- read.table("edges.txt",sep="\t",header = T) # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=26326c233918419f8dc80e8af984cdae (Download Cytoscape Data, networkedges_selfloop)
motifs <- read.csv("MS2LDA_prob0.01_overlap0.3.csv") # downloaded from http://ms2lda.org/

motifs$Document <- gsub('.*\\_', '', motifs$Document)

shmotifs <- c()
for (i in 1:nrow(edges)){
  nodes <- edges[i,c(1:2)]
  a <- motifs$Motif[motifs$Document %in% nodes[1]]
  b <- motifs$Motif[motifs$Document %in% nodes[2]]
  shmotifs <- c(shmotifs,paste(sort(intersect(a,b)),collapse = ","))
}

edges_sh <- cbind(edges,shmotifs)

compindex <- unique(edges_sh$ComponentIndex)
compindex <- compindex[-(which(compindex == -1))]

edges_new <- edges_sh
edges_new$topX <- rep(0,nrow(edges))
edges_new$topXinNode <- rep(0,nrow(edges))

for (i in 1:length(compindex)){
  pre <- c()
  clmotifs <- edges_new[which(edges_new$ComponentIndex==compindex[i]),]
  mpool <- count(sort(unlist(strsplit(as.character(clmotifs$shmotifs),split = ","))))
  topX <- as.character(mpool[order(mpool$freq,decreasing=T),][1:5,1])
  edges_new$topX[which(edges_new$ComponentIndex==compindex[i])] <- rep(paste(topX,collapse = ","),nrow(clmotifs))
  for (j in 1:nrow(clmotifs)){
    pre <- c(pre,paste(unlist(strsplit(as.character(clmotifs$shmotifs[j]),","))[which(unlist(strsplit(as.character(clmotifs$shmotifs[j]),",")) %in% topX)],collapse = ","))
  }
  edges_new$topXinNode[which(edges_new$ComponentIndex==compindex[i])] <- pre
}

# create multiple edges

edges$interact <- rep("cosine",nrow(edges))
edges_new2 <- edges[,c("CLUSTERID1","interact","CLUSTERID2","DeltaMZ","MEH","Cosine","OtherScore","ComponentIndex")]

for (i in 1:nrow(edges)){
  nodes <- edges[i,c(1:2)]
  a <- motifs$Motif[motifs$Document %in% nodes[1]]
  b <- motifs$Motif[motifs$Document %in% nodes[2]]
  c <- paste(sort(intersect(a,b)))
  
  if (length(c) != 0){
    for (j in 1:length(c)){
      edges_new2 <- rbind(edges_new2,c("CLUSTERID1"=edges[i,1],"interact"=c[j],edges[i,2:(ncol(edges)-1)]))
    }
  }
  
}

edges_final <- merge(edges_new2[edges_new2$interact=="cosine",],edges_new)
edges_final2 <- rbind.fill(edges_final,edges_new2[edges_new2$interact!="cosine",])
edges_final2 <- edges_final2[,c("CLUSTERID1", "interact", "CLUSTERID2", "DeltaMZ", "MEH", "Cosine", 
                                "OtherScore", "ComponentIndex", "shmotifs", "topX", "topXinNode")]

write.table(edges_final2,"Mass2Motifs_Edges.tsv",quote=F,row.names = F,sep="\t")

##############################################
#   map Mass2Motifs on nodes                 #
##############################################

motifs <- read.csv("MS2LDA_prob0.01_overlap0.3.csv") # downloaded from http://ms2lda.org/basicviz/short_summary/390/

motifs$Document <- gsub('.*\\_', '', motifs$Document)

motifs[-1] = apply(motifs[-1],2,as.character)
motifs_cytoscape <- aggregate(motifs[-1],by=list(motifs$Document),c)

motifs_cytoscape$Motif <-  unlist(lapply(motifs_cytoscape$Motif,paste,collapse=","))
motifs_cytoscape$Probability <-  unlist(lapply(motifs_cytoscape$Probability,paste,collapse=","))
motifs_cytoscape$Overlap.Score <-  unlist(lapply(motifs_cytoscape$Overlap.Score,paste,collapse=","))
motifs_cytoscape$Precursor.Mass <-  unlist(lapply(motifs_cytoscape$Precursor.Mass,paste,collapse=","))
motifs_cytoscape$Retention.Time <-  unlist(lapply(motifs_cytoscape$Retention.Time,paste,collapse=","))
motifs_cytoscape$Document.Annotation <-  unlist(lapply(motifs_cytoscape$Document.Annotation,paste,collapse=","))

splitmot = unique(unlist(strsplit(motifs_cytoscape$Motif, ",")))

mat <- matrix("0.00",nrow(motifs_cytoscape),length(splitmot))
colnames(mat) <- splitmot

for (i in 1:nrow(motifs_cytoscape)){
  w <- match(unlist(strsplit(motifs_cytoscape$Motif[i],",")), colnames(mat))
  mat[i,w] <- unlist(strsplit(motifs_cytoscape$Overlap.Score[i],","))
}

mat <- cbind(motifs_cytoscape,mat)
colnames(mat)[1] <- "shared.name"

write.table(mat,"Mass2Motifs_Nodes.tsv",quote=F,row.names = F,sep=";")

# import Mass2Motifs_Edges.tsv as network into Cytoscape version 3.4.0 [Shannon et al., 2003] with "CLUSTERID1" as source node, 
# "interact" as interaction type and "CLUSTERID2" as target node.

# import Mass2Motifs_Nodes.tsv as table into Cytoscape version 3.4.0 [Shannon et al., 2003]. Select semicolon as delimiter (Advanced Options) 
# and import non-numeric data columns as list of strings (data type), and comma as list delimiter.
