######################################################################################
#                                                                                    #
#      Phylogenetic generalized least squares regression analysis (PGLS)             #
#                                                                                    #
######################################################################################

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

##########################################################################
# Retrieve bioactivity data and number of nodes per chemical subclass    #
##########################################################################

## read biological activity data
biodata <- read.csv("PGLS/Network_Bioactivity.csv",sep=";")

var <- biodata[,c(2:3)]
var <- var[-which(is.na(var$Species)),]
var$Species <- as.character(var$Species)
var$Species[which(var$Species =="characias")] <- "characias.wulfenii"
var <- var[-which(var$Species =="cyparissias")[2],]

## read mass spectral molecular network data
net <- read.csv("Cytoscape_SummaryTable.csv") # Table exported from Cytoscape, with all network attributes loaded

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

classy <- net[,which(colnames(net) %in% c("componentindex","Subclass_Interpretation","Diterpenoid_Subclass_Interpretation","EuphorbiaDiterpenoid_Subclass_Interpretation"))]
classy$Diterpenoid_Subclass_Interpretation[classy$Diterpenoid_Subclass_Interpretation=="Euphorbia diterpenoid"] <- "Euphorbia diterpenoids"
classy$Diterpenoid_Subclass_Interpretation[classy$Diterpenoid_Subclass_Interpretation=="Regular diterpenoid"] <- "Regular diterpenoids"

# Chemical subclasses
cl <- sort(as.character(unique(classy$Subclass_Interpretation)))
cl <- cl[-which(cl=="-")]
cl <- sort(cl)
cl <- cl[c(2:10,1,11:12)]

subcl <- list()

for (i in 1:length(cl)){
  subcl[[i]] <- classy[which(classy$Subclass_Interpretation==cl[i]),c("componentindex","Subclass_Interpretation")]
}
names(subcl) <- cl

# Diterpenoid subclasses (Euphorbia diterpenoid, Euphorbia diterpenoids, Regular diterpenoid, Regular diterpenoids)
dit_cl <- sort(as.character(unique(classy$Diterpenoid_Subclass_Interpretation)))
dit_cl <- dit_cl[-which(dit_cl=="")]

subcldit <- list()

for (i in 1:length(dit_cl)){
  subcldit[[i]] <- classy[which(classy$Diterpenoid_Subclass_Interpretation==dit_cl[i]),c("componentindex","Diterpenoid_Subclass_Interpretation")]
}
names(subcldit) <- dit_cl

# Euphorbia diterpenoid subclasses (Phorboids, Non-phorboids)
edit_cl <- sort(as.character(unique(classy$EuphorbiaDiterpenoid_Subclass_Interpretation)))
edit_cl <- edit_cl[-which(edit_cl=="")]

subcledit <- list()

for (i in 1:length(edit_cl)){
  subcledit[[i]] <- classy[which(classy$EuphorbiaDiterpenoid_Subclass_Interpretation==edit_cl[i]),c("componentindex","EuphorbiaDiterpenoid_Subclass_Interpretation")]
}
names(subcledit) <- edit_cl

mat <- matrix(0.1,nrow(var),length(subcl)+length(c(dit_cl,edit_cl)))
colnames(mat) <- c(cl,c(dit_cl,edit_cl))
var <- cbind(var,mat)

for (i in 1:length(subcl)){
  
  ci <- as.numeric(as.character(unique(subcl[[i]]$componentindex)))
  net_red <- net[which(net[,which(colnames(net)=="componentindex")] %in% ci),]
  net_red <- net_red[,which(colnames(net_red) %in% c("componentindex",samp))]
  
  colnames(net_red)[colnames(net_red)=="GG1P7_41characias.MS1_GG1_01_28816"] <- "characias.wulfenii"
  colnames(net_red)[colnames(net_red)=="BA2P7_12acanthot_BA2_01_28605"] <- "acanthothamnos"
  
  for (j in 1:length(var$Species)){
    x <- grep(var$Species[j],colnames(net_red))
    n <- length(which(net_red[,x]!=0))
    var[j,which(colnames(var)==names(subcl)[i])]<- n
  }
}

for (i in 1:length(subcldit)){
  
  ci <- as.numeric(as.character(unique(subcldit[[i]]$componentindex)))
  net_red <- net[which(net[,which(colnames(net)=="componentindex")] %in% ci),]
  net_red <- net_red[,which(colnames(net_red) %in% c("componentindex",samp))]
  
  colnames(net_red)[colnames(net_red)=="GG1P7_41characias.MS1_GG1_01_28816"] <- "characias.wulfenii"
  colnames(net_red)[colnames(net_red)=="BA2P7_12acanthot_BA2_01_28605"] <- "acanthothamnos"
  
  for (j in 1:length(var$Species)){
    x <- grep(var$Species[j],colnames(net_red))
    n <- length(which(net_red[,x]!=0))
    var[j,which(colnames(var)==names(subcldit)[i])]<- n
  }
}

for (i in 1:length(subcledit)){
  
  ci <- as.numeric(as.character(unique(subcledit[[i]]$componentindex)))
  net_red <- net[which(net[,which(colnames(net)=="componentindex")] %in% ci),]
  net_red <- net_red[,which(colnames(net_red) %in% c("componentindex",samp))]
  
  colnames(net_red)[colnames(net_red)=="GG1P7_41characias.MS1_GG1_01_28816"] <- "characias.wulfenii"
  colnames(net_red)[colnames(net_red)=="BA2P7_12acanthot_BA2_01_28605"] <- "acanthothamnos"
  
  for (j in 1:length(var$Species)){
    x <- grep(var$Species[j],colnames(net_red))
    n <- length(which(net_red[,x]!=0))
    var[j,which(colnames(var)==names(subcledit)[i])]<- n
  }
}

var[,2:ncol(var)] <- sapply(var[,2:ncol(var)], as.character)
var[,2:ncol(var)] <- sapply(var[,2:ncol(var)], as.numeric)
rownames(var) <- var$Species
colnames(var) <- gsub(" ","",colnames(var))
colnames(var) <- gsub(";","",colnames(var))

# write dataframe to file 
write.csv(var,"PGLS/ChemicalSubclasses.csv",quote=F,row.names = F)

##########################
# Perform PGLS analysis  #
##########################

chemsub <- read.csv("PGLS/ChemicalSubclasses.csv",check.names = F)
colnames(chemsub)[colnames(chemsub)=="Terpenoids(di-,ortriterpenoids)"] <- "Terpenoids"
colnames(chemsub)[colnames(chemsub)=="Non-phorboids"] <- "Nonphorboids"

comp.data<-comparative.data(btree, chemsub, names.col="Species", vcv=T, vcv.dim=3, warn.dropped=TRUE)

# select only chemical sublcasses with at least non-zero values for at least 20 species

sort(names(which(colSums(chemsub != 0)>20)))

fit2_Benzoicacidsandderivatives<-pgls(Benzoicacidsandderivatives~Bioactive1, data=comp.data)
fit2_Cholestanesteroids<-pgls(Cholestanesteroids~Bioactive1, data=comp.data)
fit2_Diterpenoids<-pgls(Diterpenoids~Bioactive1, data=comp.data)
fit2_Glycosphingolipids<-pgls(Glycosphingolipids~Bioactive1, data=comp.data)
fit2_Glycosylglycerols<-pgls(Glycosylglycerols~Bioactive1, data=comp.data)
fit2_Stigmastanesandderivatives<-pgls(Stigmastanesandderivatives~Bioactive1, data=comp.data)
fit2_Triterpenoids<-pgls(Triterpenoids~Bioactive1, data=comp.data)
fit2_EuphorbiaDiterpenoids<-pgls(Euphorbiaditerpenoids~Bioactive1, data=comp.data)
fit2_Regularditerpenoids<-pgls(Regularditerpenoids~Bioactive1, data=comp.data)

summary(fit2_Benzoicacidsandderivatives)
summary(fit2_Cholestanesteroids)
summary(fit2_Diterpenoids) #
summary(fit2_Glycosphingolipids)
summary(fit2_Glycosylglycerols) #
summary(fit2_Stigmastanesandderivatives)
summary(fit2_Triterpenoids) 
summary(fit2_EuphorbiaDiterpenoids) #
summary(fit2_Regularditerpenoids)

chemsub <- chemsub[,which(colnames(chemsub) %in% c("Species","Bioactive1","Diterpenoids","Glycosylglycerols","Euphorbiaditerpenoids"))]

pltree <- chronopl(btree, lambda = 0, age.min = 1)
for (i in 3:ncol(chemsub)){
  dat <- matrix(0,nrow(chemsub),2)
  rownames(dat) <- chemsub$Species
  dat[,1] <- as.numeric(as.character(chemsub$Bioactive1))
  dat[,2] <- chemsub[,i]
  dat <- log(dat+1)
  dat[dat==Inf|dat== -Inf] <- 0
  dat[,1] <- dat[,1]*(-1)
  
  par(bg = "transparent") # makes background transparent
  pdf(paste(colnames(chemsub)[i], ".pdf",sep=""),width=12,height=10, bg="transparent")
  plotTree.barplot(pltree,dat,args.barplot=list(beside=TRUE,xlim=c(min(dat),max(dat)),xlab="",
                                                legend.text=TRUE,space=c(0,0.6),args.legend=list(x=max(dat),y=17,border = c("blue","red"),bty="n"),col=c("red","blue"),border=c("red","blue")),args.plotTree=list(fsize=1,color="black"))
  dev.off()
}

