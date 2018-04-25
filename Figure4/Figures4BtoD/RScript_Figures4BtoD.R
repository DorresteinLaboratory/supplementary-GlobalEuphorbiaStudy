########################################################
#                                                      #
#    R script used for creating Figures 4B to D        #
#                                                      #
########################################################

# Number of putatively annotated Euphorbia diterpenoids per species analysed
# Number of TNF-a modulating fractions per species analysed and 
# Euphorbia phylogenetic tree

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

##########################
# Plot  data             #
##########################

chemsub <- read.csv("PGLS/ChemicalSubclasses.csv",check.names = F)
colnames(chemsub)[colnames(chemsub)=="Terpenoids(di-,ortriterpenoids)"] <- "Terpenoids"
colnames(chemsub)[colnames(chemsub)=="Non-phorboids"] <- "Nonphorboids"

chemsub <- chemsub[,which(colnames(chemsub) %in% c("Species","Bioactive1","Euphorbiaditerpenoids"))]

pltree <- chronopl(btree, lambda = 0, age.min = 1)

dat <- matrix(0,nrow(chemsub),2)
rownames(dat) <- chemsub$Species
dat[,1] <- as.numeric(as.character(chemsub$Bioactive1))
dat[,2] <- chemsub[,3]
dat <- log(dat+1)
dat[dat==Inf|dat== -Inf] <- 0
dat[,1] <- dat[,1]*(-1)
  
par(bg = "transparent") # makes background transparent
pdf("Figures4BtoD",width=12,height=10, bg="transparent")
plotTree.barplot(pltree,dat,args.barplot=list(beside=TRUE,xlim=c(min(dat),max(dat)),xlab="",
                                                legend.text=TRUE,space=c(0,0.6),args.legend=list(x=max(dat),y=17,border = c("blue","red"),bty="n"),col=c("red","blue"),border=c("red","blue")),args.plotTree=list(fsize=1,color="black"))
dev.off()

