#################################################
#                                               #
#    R script used for creating Figure 1        #
#                                               #
#################################################

##############
#  Fig. 1A   #
##############

# Distribution of specialized metabolite classes on the Euphorbia phylogenetic tree

## load libraries
library(phytools)
library(caper)
library(ape)
library(geiger)
library(nlme)

## read tree
tree <- read.nexus("BayesianPhylogeneticAnalysis/EuphorbiaMBPF15_HornSet.con.tre")  

outgroup <- c("Nealchornea",
              "Anthostema",
              "Bonania",
              "Calycopeplus",
              "Colliguaja",
              "Dichostemma",
              "Gymnanthes",
              "Homalanthus",
              "Hura",
              "Mabea",
              "Maprounea",
              "Microstachys",
              "Neoguillauminia",
              "Senefelderopsis",
              "Stillingia")

## root tree
tree <- root(tree,outgroup = outgroup,resolve.root = T)
tree <- drop.tip(tree,as.character(tree$tip.label[which(tree$tip.label %in% outgroup)]))

plotTree(tree,fsize=0.7,node.numbers=T)

## add tips to tree 
# sarcoceras Alectoroctonum with cotinifolia N58
# hyssopifolia Anisophyllum with hirta N57
# ophthalmica Anisophyllum with hirta N57
# prostrata Anisophyllum with hirta N57
# thymifolia Anisophyllum with hirta N57
# [Yang et al., 2012]

btree<-bind.tip(tree,"sarcoceras",where=which(tree$tip.label=="cotinifolia"),position= 0.5*tree$edge.length[which(tree$edge[,2]==which(tree$tip.label=="cotinifolia"))],edge.length=tree$edge.length[which(tree$edge[,2]==which(tree$tip.label=="cotinifolia"))])
btree<-bind.tip(btree,"hyssopifolia",where=which(btree$tip.label=="hirta"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))])
btree<-bind.tip(btree,"ophthalmica",where=which(btree$tip.label=="hyssopifolia"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))])
btree<-bind.tip(btree,"prostrata",where=which(btree$tip.label=="hyssopifolia"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))])
btree<-bind.tip(btree,"thymifolia",where=which(btree$tip.label=="hirta"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))])

btree$tip.label[btree$tip.label=="peplus.2"] <- "peplus"
btree$tip.label[btree$tip.label=="horrida.1"] <- "horrida"
btree$tip.label[btree$tip.label=="davidii"] <- "dentata"

plotTree(btree,fsize=0.7,node.numbers=T)

## read mass spectral molecular network data
net <- read.csv("Cytoscape_SummaryTable.csv") # Table exported from Cytoscape, with all network attributes loaded
## define sample names
samp <- c("BA2P7_12acanthot_BA2_01_28605",                              
          "BC11P7_211aeruginosa_BC11_01_28719",                          
          "BC1P7_21helioscopia_BC1_01_28669",                           
          "BC2P7_22lagascae_BC2_01_28674",                               
          "BC7P7_27ornithopus_BC7_01_28699",                            
          "BC8P7_28balsamifera_BC8_01_28704",     
          "GA10P7_110lactea.MS1_GA10_01_28647",                         
          "GA11P7_111cylindrifolia.MS1_GA11_01_28652",                   
          "GA12P7_112sipolisii.MS1_GA12_01_28657",                      
          "GA1P7_11peplus.MS1_GA1_01_28602",                             
          "GA3P7_13cyparissias.MS1_GA3_01_28612",                       
          "GA4P7_14platyclada.MS1_GA4_01_28617",                         
          "GA5P7_15hyssopifolia.MS1_GA5_01_28622",                      
          "GA6P7_16ophthalmica.MS1_GA6_01_28627",                        
          "GA7P7_17mammillaris.MS1_GA7_01_28632",                       
          "GA8P7_18obesa.MS1_GA8_01_28637",                              
          "GA9P7_19grandicornis.MS1_GA9_01_28642",                      
          "GC10P7_21weberbaueri.MS1_GC10_01_28716",                      
          "GC12P7_212abdelkuri.MS1.RERUN_GC12_01_28881",                
          "GC3P7_23amygdaloides.MS1_GC3_01_28681",                       
          "GC4P7_24cotinifolia.MS1_GC4_01_28686",                       
          "GC5P7_25hirta.MS1_GC5_01_28691",                              
          "GC6P7_26prostrata.MS1_GC6_01_28696",                         
          "GC9P7_29ingens.MS1_GC9_01_28711",                             
          "GE10P7_31stenoclada.MS1_GE10_01_28789",                      
          "GE11P7_311ammak.MS1.RERUN_GE11_01_28883",                     
          "GE12P7_312fiherenensis.MS1_GE12_01_28811",                   
          "GE1P7_31nicaeensis.MS1_GE1_01_28744",                         
          "GE2P7_32lathyris.MS1_GE2_01_28749",                          
          "GE3P7_33myrsinites.MS1_GE3_01_28754",                         
          "GE4P7_34dentata.MS1_GE4_01_28759",                           
          "GE5P7_35sarcoceras.MS1_GE5_01_28764",                         
          "GE6P7_36graminea.MS1_GE6_01_28769",                          
          "GE7P7_37globosa.MS1_GE7_01_28774",                            
          "GE8P7_38horrida.MS1_GE8_01_28779",                           
          "GE9P7_39milii.MS1_GE9_01_28784",                              
          "GG1P7_41characias.MS1_GG1_01_28816",                         
          "GG2P7_42segetalis.MS1_GG2_01_28821",                          
          #"GG3P7_43cyparissias.MS1_GG3_01_28826"                       
          "GG4P7_44thymifolia.MS1_GG4_01_28831",                         
          "GG5P7_45bubalina.MS1_GG5_01_28836",                          
          "GG6P7_46jansenvillensis.MS1_GG6_01_28841",                    
          "GG7P7_47neriifolia.MS1_GG7_01_28846",                        
          "GG8P7_48alluaudii.MS1_GG8_01_28851")

net_samp <- net[,which(colnames(net) %in% c("componentindex",samp))]
colnames(net_samp)[colnames(net_samp)=="GG1P7_41characias.MS1_GG1_01_28816"] <- "characias.wulfenii"
colnames(net_samp)[colnames(net_samp)=="BA2P7_12acanthot_BA2_01_28605"] <- "acanthothamnos"

classy <- net[,which(colnames(net) %in% c("componentindex","Subclass_Interpretation"))]

cl <- sort(as.character(unique(classy$Subclass_Interpretation)))
#cl[which(cl=="\"Terpenoids (di-, or triterpenoids)\"")] <- "Terpenoids (di-, or triterpenoids)"
cl <- cl[-which(cl=="-")]
cl <- sort(cl)
cl <- cl[c(2:10,1,11:12)]

mat <- c()
for (i in 1:length(cl)){
  comp <- classy[which(classy$Subclass_Interpretation==cl[i]),]
  
  #diter_red <- diter[-which(diter$CF_directparent %in% noneuphdit),]
  
  ci <- as.numeric(as.character(unique(comp$componentindex)))
  if (length(which(is.na(ci)==T))!=0){
    ci <- ci[-which(is.na(ci)==T)]
  }
  
  nodes <- net_samp[which(net_samp[,which(colnames(net_samp)=="componentindex")] %in% ci),]
  
  nums <- c()
  for (j in 1:length(btree$tip.label)){
    x <- grep(btree$tip.label[j],colnames(nodes))
    n <- length(which(nodes[,x]!=0))
    nums <- c(nums,n)
  }
  mat <- cbind(mat,nums)
}

colnames(mat) <- cl
rownames(mat) <- btree$tip.label
mat[mat==0] <- NA

colors<-colorRampPalette(colors=c("#C4DBC4","darkgreen"))(6)

pdf("Figure1A.pdf",width=7,height=10, bg="transparent")
phylo.heatmap(btree,log(mat),fsize=0.8,colors=colors,standardize = F)
dev.off()

##############
#  Fig. 1B   #
##############

# Molecular features representing individual mass spectral molecular network nodes shared across 
# species of Euphorbia subgeneric clades

## load libraries
library(VennDiagram)
library(RColorBrewer)
library(venneuler)

## read mass spectral molecular network data
net <- read.csv("Cytoscape_SummaryTable.csv") # Table exported from Cytoscape, with all network attributes loaded

cols <- c("shared.name","FTIC.SUBG_Athymalus","FTIC.SUBG_Chamaesyce","FTIC.SUBG_Esula","FTIC.SUBG_Euphorbia")
net <- net[,which(colnames(net)%in%cols)]

## create Venn Diagram

VENN.LIST <- list()

for (i in 1:4){
  VENN.LIST[[i]] <- as.character(net$shared.name)[which(net[,i]!=0)]
}

names(VENN.LIST) <- c("Athymalus","Chamaesyce","Esula","Euphorbia")

fillcols <- brewer.pal(n = 4, name = "Set1")

venn.plot <- venn.diagram(VENN.LIST , NULL, fill=fillcols, alpha=c(0.7,0.7,0.7,0.7), 
                          cex = 2, cat.fontface=4, category.names="", 
                          main="",lty =1,cat.fontfamily = "sans")

# plot Venn Diagram (accurate)
grid.draw(venn.plot)

# plot Venn Diagram with proportional size (approximate)
v <- venneuler(c(A=2018, B=540, C=3584, D=2095, "A&B"=17, "A&C"=299,"A&D"=249,"A&C&D"=427,"A&B&D"=23,"A&B&C"=23,"A&B&C&D"=207,
                 "C&D"=312, "B&C&D"=73,"B&C"=58,"B&D"=37))

pdf("VennDiagram_Proportional.pdf",width=10,height=12, bg="transparent") 
plot(v,col=fillcols)
dev.off()

##############
#  Fig. 1C   #
##############

# Chemical similarity among Euphorbia subgeneric clades assessed using the chemical structural compositional
# similarity [Sedio et al., 2017]

## read mass spectral molecular network data
net <- read.csv("Cytoscape_SummaryTable.csv") # Table exported from Cytoscape, with all network attributes loaded

samples_athym <- c("BC7P7_27ornithopus_BC7_01_28699","BC8P7_28balsamifera_BC8_01_28704",
                   "GA7P7_17mammillaris.MS1_GA7_01_28632","GA8P7_18obesa.MS1_GA8_01_28637","GE7P7_37globosa.MS1_GE7_01_28774",
                   "GE8P7_38horrida.MS1_GE8_01_28779","GG5P7_45bubalina.MS1_GG5_01_28836","GG6P7_46jansenvillensis.MS1_GG6_01_28841")

samples_cham <- c("GA4P7_14platyclada.MS1_GA4_01_28617","GA5P7_15hyssopifolia.MS1_GA5_01_28622","GA6P7_16ophthalmica.MS1_GA6_01_28627",
                  "GC4P7_24cotinifolia.MS1_GC4_01_28686","GC5P7_25hirta.MS1_GC5_01_28691","GC6P7_26prostrata.MS1_GC6_01_28696",
                  "GE4P7_34dentata.MS1_GE4_01_28759","GE5P7_35sarcoceras.MS1_GE5_01_28764","GE6P7_36graminea.MS1_GE6_01_28769",
                  "GG4P7_44thymifolia.MS1_GG4_01_28831")

samples_esula <- c("GA1P7_11peplus.MS1_GA1_01_28602","GA3P7_13cyparissias.MS1_GA3_01_28612","GC3P7_23amygdaloides.MS1_GC3_01_28681",
                   "GE1P7_31nicaeensis.MS1_GE1_01_28744","GE2P7_32lathyris.MS1_GE2_01_28749","GE3P7_33myrsinites.MS1_GE3_01_28754",
                   "GG2P7_42segetalis.MS1_GG2_01_28821","GG1P7_41characias.MS1_GG1_01_28816",
                   "BA2P7_12acanthot_BA2_01_28605","BC2P7_22lagascae_BC2_01_28674","BC1P7_21helioscopia_BC1_01_28669")

samples_euph <- c("BC11P7_211aeruginosa_BC11_01_28719","GA10P7_110lactea.MS1_GA10_01_28647","GA11P7_111cylindrifolia.MS1_GA11_01_28652",
                  "GA12P7_112sipolisii.MS1_GA12_01_28657","GA9P7_19grandicornis.MS1_GA9_01_28642","GC10P7_21weberbaueri.MS1_GC10_01_28716",
                  "GC12P7_212abdelkuri.MS1.RERUN_GC12_01_28881","GC9P7_29ingens.MS1_GC9_01_28711","GE10P7_31stenoclada.MS1_GE10_01_28789",
                  "GE11P7_311ammak.MS1.RERUN_GE11_01_28883","GE12P7_312fiherenensis.MS1_GE12_01_28811",
                  "GG7P7_47neriifolia.MS1_GG7_01_28846","GG8P7_48alluaudii.MS1_GG8_01_28851","GE9P7_39milii.MS1_GE9_01_28784")

features <- net[,which(colnames(net) %in% c(samples_athym,samples_cham,samples_esula,samples_euph))]
rownames(features) <- net$shared.name

## calculate chemical structural compositional similarity per Euphorbia subgeneric clade

## change according to subgenus
tab <- features[,which(colnames(features) %in% c("shared.name",samples_athym))]
## change according to subgenus
tab <- t(apply(t(tab), 1, function(x) x/sum(as.numeric(x))))

samp_pairs <- combn((rownames(tab)),2)
cos <- read.delim('Mass2Motifs_2_MolecularNetwork/edges.txt') # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=26326c233918419f8dc80e8af984cdae (Download Cytoscape Data, networkedges_selfloop)s
cos[cos[,5] < 0.6,5] <- 0

cscs <- c()

for (i in 1:ncol(samp_pairs)){
  
  ab <- matrix(tab[which(rownames(tab)==samp_pairs[1,i]),], ncol=1) %*% matrix(tab[which(rownames(tab)==samp_pairs[2,i]),], nrow=1)
  aa <- matrix(tab[which(rownames(tab)==samp_pairs[1,i]),], ncol=1) %*% matrix(tab[which(rownames(tab)==samp_pairs[1,i]),], nrow=1)
  bb <- matrix(tab[which(rownames(tab)==samp_pairs[2,i]),], ncol=1) %*% matrix(tab[which(rownames(tab)==samp_pairs[2,i]),], nrow=1)
  
  dcos <- matrix(0, nrow=ncol(tab), ncol=ncol(tab))
  
  dcos[cbind(match(cos[,1], colnames(tab)), match(cos[,2], colnames(tab)))] <- cos[,5]
  dcos[cbind(match(cos[,2], colnames(tab)), match(cos[,1], colnames(tab)))] <- cos[,5]
  diag(dcos) <- 1
  
  cosXaa = dcos * aa
  cosXbb = dcos * bb
  cosXab = dcos * ab
  
  score <- sum(cosXab)/max(c(sum(cosXaa), sum(cosXbb)))
  
  cscs <- c(cscs,score)
}

cscs_pairs <- rbind(samp_pairs,cscs)
cscs_pairs <- t(cscs_pairs)
colnames(cscs_pairs)[1:2] <- c("sample1","sample2")

write.csv(cscs_pairs,"cscs_athymalus.txt",quote=F,row.names = F)
write.csv(cscs_pairs,"cscs_chamaesyce.txt",quote=F,row.names = F)
write.csv(cscs_pairs,"cscs_esula.txt",quote=F,row.names = F)
write.csv(cscs_pairs,"cscs_euphorbia.txt",quote=F,row.names = F)

################ plot results

## load libraries
library(ggplot2)
library(plyr)
library(multcompView)

cscs_athym_tab <- read.csv("cscs_athymalus.txt")
cscs_athym <- cscs_athym_tab$cscs

cscs_cham_tab <- read.csv("cscs_chamaesyce.txt")
cscs_cham <- cscs_cham_tab$cscs

cscs_esula_tab <- read.csv("cscs_esula.txt")
cscs_esula <- cscs_esula_tab$cscs

cscs_euph_tab <- read.csv("cscs_euphorbia.txt")
cscs_euph <- cscs_euph_tab$cscs

df <- c(cscs_athym, cscs_cham, cscs_esula,cscs_euph)
df_n <- c(rep("Athymalus",length(cscs_athym)),rep("Chamaesyce",length(cscs_cham)),rep("Esula",length(cscs_esula)),rep("Euphorbia",length(cscs_euph)))

df_t <- cbind(df,df_n)
df_t <- as.data.frame(df_t)
df_t$df <- as.numeric(as.character(df_t$df))
cols <- c('#e41a1c','#377eb8', '#4daf4a','#984ea3')
df_cols <- cols[as.numeric(df_t$df_n)]

p<-ggplot(df_t, aes(x=df_n, y=df, color=df_n)) +
  geom_boxplot(color='black',fill=cols) +
  geom_point(size=1,color='black')   +
  xlab("")+
  ylab("Chemical structural compositional similarity, CSCS") +
  scale_y_continuous(limits = c(0, 1))
p

# Compute and plot Tukey Honest Significant Differences

m <- aov(df ~ df_n,data=df_t)
tHSD <- TukeyHSD(m, ordered = FALSE, conf.level = 0.95)

# https://stackoverflow.com/questions/44712185/tukeys-post-hoc-on-ggplot-boxplot?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
####################################################################### Tukey lable plotting function
generate_label_df <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(df_t, flev, function (x) max(fivenum(x$df)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}
####################################################################### Tukey lable plotting function

pdf("Figure1C.pdf",width=4,height=7, bg="transparent") 
p <- ggplot(df_t,aes(x=df_n,y=df))+ 
  geom_boxplot(color='black',fill=cols) +
  geom_text(data = generate_label_df(tHSD, 'df_n'), aes(x = plot.labels, y = V1, label = labels))+
  xlab("")+
  ylab("Chemical structural compositional similarity (CSCS)") +
  geom_point(size=1.5,color='black')   +
  scale_y_continuous(limits = c(0, 1.2)) +
  theme_light()
p
dev.off()

##############
#  Fig. 1D   #
##############

# Euphorbia chemogram (hierarchical cluster analysis of the pair-wise chemical structural compositional dissimilarities)
# and phylogenetic tree

## load libraries
library(phytools)
library(caper)
library(ape)
library(geiger)
library(nlme)
library(dendextend)

## read tree
tree <- read.nexus("BayesianPhylogeneticAnalysis/EuphorbiaMBPF15_HornSet.con.tre")  

outgroup <- c("Nealchornea",
              "Anthostema",
              "Bonania",
              "Calycopeplus",
              "Colliguaja",
              "Dichostemma",
              "Gymnanthes",
              "Homalanthus",
              "Hura",
              "Mabea",
              "Maprounea",
              "Microstachys",
              "Neoguillauminia",
              "Senefelderopsis",
              "Stillingia")

## root tree
tree <- root(tree,outgroup = outgroup,resolve.root = T)
tree <- drop.tip(tree,as.character(tree$tip.label[which(tree$tip.label %in% outgroup)]))

plotTree(tree,fsize=0.7,node.numbers=T)

## add tips to tree 
# sarcoceras Alectoroctonum with cotinifolia N58
# hyssopifolia Anisophyllum with hirta N57
# ophthalmica Anisophyllum with hirta N57
# prostrata Anisophyllum with hirta N57
# thymifolia Anisophyllum with hirta N57
# [Yang et al., 2012]

btree<-bind.tip(tree,"sarcoceras",where=which(tree$tip.label=="cotinifolia"),position= 0.5*tree$edge.length[which(tree$edge[,2]==which(tree$tip.label=="cotinifolia"))],edge.length=tree$edge.length[which(tree$edge[,2]==which(tree$tip.label=="cotinifolia"))])
btree<-bind.tip(btree,"hyssopifolia",where=which(btree$tip.label=="hirta"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))])
btree<-bind.tip(btree,"ophthalmica",where=which(btree$tip.label=="hyssopifolia"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))])
btree<-bind.tip(btree,"prostrata",where=which(btree$tip.label=="hyssopifolia"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hyssopifolia"))])
btree<-bind.tip(btree,"thymifolia",where=which(btree$tip.label=="hirta"),position= 0.5*btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))],edge.length=btree$edge.length[which(btree$edge[,2]==which(btree$tip.label=="hirta"))])

btree$tip.label[btree$tip.label=="peplus.2"] <- "peplus"
btree$tip.label[btree$tip.label=="horrida.1"] <- "horrida"
btree$tip.label[btree$tip.label=="davidii"] <- "dentata"

plotTree(btree,fsize=0.7,node.numbers=T)

###################### create chemogram

## read mass spectral molecular network data
net <- read.csv("Cytoscape_SummaryTable.csv") # Table exported from Cytoscape, with all network attributes loaded

samples <- c("BC7P7_27ornithopus_BC7_01_28699","BC8P7_28balsamifera_BC8_01_28704",
             "GA7P7_17mammillaris.MS1_GA7_01_28632","GA8P7_18obesa.MS1_GA8_01_28637","GE7P7_37globosa.MS1_GE7_01_28774",
             "GE8P7_38horrida.MS1_GE8_01_28779","GG5P7_45bubalina.MS1_GG5_01_28836","GG6P7_46jansenvillensis.MS1_GG6_01_28841",
             "GA4P7_14platyclada.MS1_GA4_01_28617","GA5P7_15hyssopifolia.MS1_GA5_01_28622","GA6P7_16ophthalmica.MS1_GA6_01_28627",
             "GC4P7_24cotinifolia.MS1_GC4_01_28686","GC5P7_25hirta.MS1_GC5_01_28691","GC6P7_26prostrata.MS1_GC6_01_28696",
             "GE4P7_34dentata.MS1_GE4_01_28759","GE5P7_35sarcoceras.MS1_GE5_01_28764","GE6P7_36graminea.MS1_GE6_01_28769",
             "GG4P7_44thymifolia.MS1_GG4_01_28831","GA1P7_11peplus.MS1_GA1_01_28602","GA3P7_13cyparissias.MS1_GA3_01_28612","GC3P7_23amygdaloides.MS1_GC3_01_28681",
             "GE1P7_31nicaeensis.MS1_GE1_01_28744","GE2P7_32lathyris.MS1_GE2_01_28749","GE3P7_33myrsinites.MS1_GE3_01_28754",
             "GG2P7_42segetalis.MS1_GG2_01_28821","GG1P7_41characias.MS1_GG1_01_28816",
             "BA2P7_12acanthot_BA2_01_28605","BC2P7_22lagascae_BC2_01_28674","BC1P7_21helioscopia_BC1_01_28669","BC11P7_211aeruginosa_BC11_01_28719",
             "GA10P7_110lactea.MS1_GA10_01_28647","GA11P7_111cylindrifolia.MS1_GA11_01_28652",
             "GA12P7_112sipolisii.MS1_GA12_01_28657","GA9P7_19grandicornis.MS1_GA9_01_28642","GC10P7_21weberbaueri.MS1_GC10_01_28716",
             "GC12P7_212abdelkuri.MS1.RERUN_GC12_01_28881","GC9P7_29ingens.MS1_GC9_01_28711","GE10P7_31stenoclada.MS1_GE10_01_28789",
             "GE11P7_311ammak.MS1.RERUN_GE11_01_28883","GE12P7_312fiherenensis.MS1_GE12_01_28811",
             "GG7P7_47neriifolia.MS1_GG7_01_28846","GG8P7_48alluaudii.MS1_GG8_01_28851","GE9P7_39milii.MS1_GE9_01_28784")

## change according to subgenus
tab <- net[,which(colnames(net) %in% c(samples))]
rownames(tab) <- net$shared.name
## change according to subgenus
tab <- t(apply(t(tab), 1, function(x) x/sum(as.numeric(x))))

cos <- read.delim('Mass2Motifs_2_MolecularNetwork/edges.txt') # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=26326c233918419f8dc80e8af984cdae (Download Cytoscape Data, networkedges_selfloop)s

source("CSCD.R")

tStart <- Sys.time() 
mcscd <- cscd(tab, cos, norm=TRUE) 
tEnd <- Sys.time() 
tEnd-tStart # Time difference of 18.96738 mins

pos <- unlist(sapply(btree$tip.label,grep,rownames(mcscd)))
rownames(mcscd)[pos] <- names(pos)
rownames(mcscd)[which(rownames(mcscd)=="GG1P7_41characias.MS1_GG1_01_28816")] <- "characias.wulfenii"
rownames(mcscd)[which(rownames(mcscd)=="BA2P7_12acanthot_BA2_01_28605")] <- "acanthothamnos"

dend <-  mcscd %>% as.dist %>%
  hclust %>% as.dendrogram 

plot(dend)

ultra_tree <- chronopl(btree,lambda = 0, age.min = 1)
dend2 <- as.dendrogram(ultra_tree)

dl <- dendlist(dend, dend2)
tanglegram(dl, sort = TRUE, common_subtrees_color_lines = F, highlight_distinct_edges  = F, highlight_branches_lwd = F)
