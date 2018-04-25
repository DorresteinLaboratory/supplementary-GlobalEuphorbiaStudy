##########################################################################################
#                                                                                        #
#   Calculate consensus chemical classes per molecular families using ClassyFire         #
#                                                                                        #
##########################################################################################

sc <- read.csv("Sirius_CSIFingerID/Cytoscape_features_Sirius.csv") 
gnps <- read.csv("gnps_node_attributes.tsv",sep="\t") # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=26326c233918419f8dc80e8af984cdae (Download Cytoscape data, clusterinfosummarygroup_attributes_withIDs_withcomponentID)
colnames(gnps)[colnames(gnps)=="cluster.index"] <- "shared.name"

sc <- merge(sc,gnps,by="shared.name",all = T)
  
ci <- unique(sc$componentindex)
ci <- ci[-which(ci==-1)] # remove single nodes

# downloaded from https://proteomics.ucsd.edu/ProteoSAFe/status.jsp?task=184a80db74334668ae1d0c0f852cb77c
naph <- read.table("node_attributes_table_H.tsv",sep="\t",quote = "",header=T)
# downloaded from https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=2cfddd3b8b1e469181a13e7d3a867a6f 
napnh4 <- read.table("node_attributes_table_NH4.tsv",sep="\t",quote = "",header=T)

naph <- naph[,c(which(colnames(naph) %in% c("cluster.index","Smiles","MetFragSMILES", "ConsensusSMILES", "FusionSMILES")))]
napnh4 <- napnh4[,c(which(colnames(napnh4) %in% c("cluster.index","Smiles","MetFragSMILES", "ConsensusSMILES", "FusionSMILES")))]

clin <- c()
sm <- c()
for (i in 1:length(ci)){
  
  sc_smiles <- sc[,c(which(colnames(sc)== "smiles"),which(colnames(sc)== "shared.name"), which(colnames(sc)== "componentindex"))]
  sc_smiles <- sc_smiles[sc_smiles$componentindex==ci[i],]
  
  siriuscsi <- unname(unlist(sapply(as.character(sc_smiles$smiles),strsplit,split=",")))
  siriuscsi <- siriuscsi[-which(siriuscsi== "no FingerID found")]
  
  NAPH <- naph[which(naph$cluster.index %in% sc_smiles$shared.name),]
  NAPNH4 <- napnh4[which(napnh4$cluster.index %in% sc_smiles$shared.name),]

  HGNPS <- unname(unlist(sapply(as.character(NAPH$Smiles),strsplit,split=",")))
  HMetFrag <- unname(unlist(sapply(as.character(NAPH$MetFragSMILES),strsplit,split=",")))
  HFusion <- unname(unlist(sapply(as.character(NAPH$FusionSMILES),strsplit,split=",")))
  HConsensus <- unname(unlist(sapply(as.character(NAPH$ConsensusSMILES),strsplit,split=",")))
  
  NH4GNPS <- unname(unlist(sapply(as.character(NAPNH4$Smiles),strsplit,split=",")))
  NH4MetFrag <- unname(unlist(sapply(as.character(NAPNH4$MetFragSMILES),strsplit,split=",")))
  NH4Fusion <- unname(unlist(sapply(as.character(NAPNH4$FusionSMILES),strsplit,split=",")))
  NH4Consensus <- unname(unlist(sapply(as.character(NAPNH4$ConsensusSMILES),strsplit,split=",")))
  
  allsmiles <- unique(c(siriuscsi,HGNPS,HMetFrag,HFusion,HConsensus,HGNPS,NH4MetFrag,NH4Fusion,NH4Consensus))
  if (length(allsmiles[which(is.na(allsmiles))]) !=0){
    allsmiles <- allsmiles[-which(is.na(allsmiles))]
  }
  
  if(length(allsmiles) != 0){
    sm <- c(sm,allsmiles)
    ident <- paste(as.character(rep(ci[i],length(allsmiles))),as.character(1:length(allsmiles)),sep = "_")
    clin <- c(clin,ident)
  }
}

smnet <- cbind(clin,sm)

write.table(smnet,"ClassyFire_InputSMILES.tsv",sep="\t",row.names = F,quote = F,col.names = F)

# ClassyFire_InputSMILES.tsv was submited to ClassyFire (http://classyfire.wishartlab.com/) [Djoumbou Feunang et al., 2016] and .json files were downloaded

###############################
# retrieve ClassyFire results #
###############################

library(rjson)

sc <- read.csv("Sirius_CSIFingerID/Cytoscape_features_Sirius.csv") 
gnps <- read.csv("gnps_node_attributes.tsv",sep="\t") # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=26326c233918419f8dc80e8af984cdae (Download Cytoscape data, clusterinfosummarygroup_attributes_withIDs_withcomponentID)
colnames(gnps)[colnames(gnps)=="cluster.index"] <- "shared.name"

sc <- merge(sc,gnps,by="shared.name",all = T)

my.JSON1 <- fromJSON(file="Network_automated_all_Part1.json") # downloaded from ClassyFire
my.JSON2 <- fromJSON(file="Network_automated_all_Part2.json") # downloaded from ClassyFire

kingdom <- list()
superclass <- list()
class <- list()
subclass <- list()
direct_parent <- list()
substituents <- list()
names_sum <- list()

for (i in 1:length(my.JSON1$entities)){
  kingdom[[i]] <- my.JSON1$entities[[i]]$kingdom$name
  superclass[[i]] <- my.JSON1$entities[[i]]$superclass$name
  class[[i]] <- my.JSON1$entities[[i]]$class$name
  subclass[[i]] <- my.JSON1$entities[[i]]$subclass$name
  direct_parent[[i]] <- my.JSON1$entities[[i]]$direct_parent$name
  substituents[[i]] <- my.JSON1$entities[[i]]$substituents
  names_sum[[i]] <- my.JSON1$entities[[i]]$identifier
  }

for (i in 1:length(my.JSON2$entities)){
  kingdom[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$kingdom$name
  superclass[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$superclass$name
  class[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$class$name
  subclass[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$subclass$name
  direct_parent[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$direct_parent$name
  substituents[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$substituents
  names_sum[[i+length(my.JSON1$entities)]] <- my.JSON2$entities[[i]]$identifier
}

kingdom <- unlist(kingdom)
superclass <- unlist(superclass)
class <- unlist(class)
wsubclass <- sapply(subclass, is.null)
subclass[wsubclass==T] <- "none"
subclass <- unlist(subclass)
direct_parent <- unlist(direct_parent)

networkgr <- sub("_[^_]+$", "", unlist(names_sum))
networks <- unique(networkgr)

clustersummary <- list()
clustersummary_scores <- list()

for (i in 1:length(networks)){
  
  w <- which(networkgr == networks[i])
  
  kingdom_num <- table(kingdom[w])/length(w)
  superclass_num <-table(superclass[w])/length(w)
  class_num <-table(class[w])/length(w)
  subclass_num <-table(subclass[w])/length(w)
  direct_parent_num <-table(direct_parent[w])/length(w)
  substituents_num <- table(unlist(substituents[w]))/length(unlist(substituents[w]))
 
  l <- list()
  l_scores <- list()
  
  l[[1]] <- names(which(kingdom_num==max(kingdom_num)))
  l[[2]] <- names(which(superclass_num==max(superclass_num)))
  l[[3]] <- names(which(class_num==max(class_num)))
  l[[4]] <- names(which(subclass_num==max(subclass_num)))
  l[[5]] <- names(which(direct_parent_num==max(direct_parent_num)))
  
  for (j in 1:length(l)){
    if(length(l[[j]]) !=1){
      l[[j]] <- paste(sort(l[[j]]),collapse = ";")
    }
    else{
      l[[j]] <- l[[j]]
    }
  }
  
  l[[6]] <- paste(sort(names(tail(sort(substituents_num),10))),collapse = ";")
  
  clustersummary[[i]] <- unlist(l)
  clustersummary_scores[[i]] <- c(round(max(kingdom_num),3),
                             round(max(superclass_num),3),
                             round(max(class_num),3),
                             round(max(subclass_num),3),
                             round(max(direct_parent_num),3),
                             paste(round(substituents_num[which(names(substituents_num) %in% sort(names(tail(sort(substituents_num),10))))],3),collapse = ";"))
  
}

df <- data.frame(matrix(unlist(clustersummary), nrow=length(clustersummary), byrow=T))
df$componentindex <- networks
colnames(df) <- c("CF_kingdom","CF_superclass","CF_class","CF_subclass","CF_directparent","CF_substituents","componentindex")
df_scores <- data.frame(matrix(unlist(clustersummary_scores), nrow=length(clustersummary_scores), byrow=T))
df_scores$componentindex <- networks
colnames(df_scores) <- c("CF_kingdom_scores","CF_superclass_scores","CF_class_scores","CF_subclass_scores","CF_directparent_scores","CF_substituents_scores","componentindex")

finalout <- merge(df,df_scores,sort=FALSE)
finalout <- merge(sc[,which(colnames(sc) %in% c("shared.name","componentindex"))],finalout,sort=FALSE,all = T)
finalout <- finalout[,c(which(colnames(finalout)=="shared.name"),which(colnames(finalout)!="shared.name"))]
colnames(finalout)[which(colnames(finalout)=="componentindex")] <- "CF_componentindex"

write.table(finalout,"ClassyFire_Output_forCytoscape.tsv",sep="\t",row.names = F,quote = F)

# import ClassyFire_Output_forCytoscape.tsv as table into Cytoscape version 3.4.0 [Shannon et al., 2003]. Import all data columns 
# as list of strings (data type), and semicolon as list delimiter.