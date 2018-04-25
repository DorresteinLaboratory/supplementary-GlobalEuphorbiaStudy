##########################################################################################
#                                                                                        #
#                 Create mapping files for visualization in 'ili                         #
#                                                                                        #
##########################################################################################

# merge metadata with fearure table
ClApp <- read.csv("MZmine/3DModels_MS2_gapfilled.csv",sep=",",header = T)
IDs <- paste(round(ClApp$row.m.z,4),round(ClApp$row.retention.time,4),paste("ID",ClApp$row.ID,sep=".."),sep="_")

ClApp <- cbind(IDs, ClApp)
ClApp <- ClApp[,-which(colnames(ClApp) %in% c("X","row.m.z","row.retention.time"))]

# remove blank samples and QCs
ClApp <- ClApp[,-grep("QC",colnames(ClApp))]
ClApp <- ClApp[,-grep("Blank",colnames(ClApp))]

ClApp <- t(ClApp)
colnames(ClApp) <- ClApp[1,]
ClApp <- ClApp[-c(1,2),]
class(ClApp) <- "numeric"

rownames(ClApp) <- gsub("X","",rownames(ClApp))
rownames(ClApp) <- gsub(".Peak.area","",rownames(ClApp))
rownames(ClApp) <- gsub("mzML","mzXML",rownames(ClApp))

# load metadata
MetaData <- read.csv("MetaData_3DMolecularCartography.csv", sep=",", dec=".", skip= 0,header = TRUE)
MetaData  <- MetaData[,-grep("X",colnames(MetaData))]
MetaData$FileName <- gsub("-",".",MetaData$FileName)

comb <- cbind(MetaData[match(rownames(ClApp),MetaData$FileName),],ClApp)
rownames(comb) <- rownames(ClApp)
ClApp <- comb

write.csv(ClApp,file="FeatureTable_MetaData_3DMolecularCartography.csv",row.names=F)

# create individual mapping files for visualization in 'ili 

#################
# Milii Model 1 #
#################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)

M1milii <- ClApp[ClApp$ModelSpecies=="M1milii",]
M1milii$PlantPartID <- gsub("roots4","roots",M1milii$PlantPartID)
colnames(M1milii)[colnames(M1milii)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EMilii_coordinates.csv", sep=";", dec=".",header = TRUE) 
coord$Sample.name <- gsub("thorns","spines",coord$Sample.name)

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M1milii$Sample.name))] # is any sample missing?

mzcoord <- merge(M1milii,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]
mzcoord2 <- mzcoord2[-which(mzcoord2$Sample.name==notin),]

write.table(mzcoord2,file="EMilii_Model1_features_MZmine.csv",row.names=F,sep=",",quote=T)

#################
# Milii Model 2 #
#################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M2milii <- ClApp[ClApp$ModelSpecies=="M2milii",]
M2milii$PlantPartID <- gsub("roots1","roots",M2milii$PlantPartID)
colnames(M2milii)[colnames(M2milii)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EMilii_coordinates.csv", sep=";", dec=".",header = TRUE) 
coord$Sample.name <- gsub("thorns","spines",coord$Sample.name)

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M2milii$Sample.name))] # is any sample missing?

mzcoord <- merge(M2milii,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="EMilii_Model2_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Lathyris Model 1 #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M1Lathyris <- ClApp[ClApp$ModelSpecies=="M1lathyris",]
colnames(M1Lathyris)[colnames(M1Lathyris)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("ELathyris_coordinates.csv", sep=";", dec=".",header = TRUE) 

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M1Lathyris$Sample.name))] # is any sample missing?

mzcoord <- merge(M1Lathyris,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="ELathyris_Model1_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Lathyris Model 2 #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M2Lathyris <- ClApp[ClApp$ModelSpecies=="M2lathyris",]
M2Lathyris$PlantPartID <- gsub("roots1","roots",M2Lathyris$PlantPartID)
colnames(M2Lathyris)[colnames(M2Lathyris)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("ELathyris_coordinates.csv", sep=";", dec=".",header = TRUE) 
coord$Sample.name <- gsub("innerroot","roots",coord$Sample.name)
coord$Sample.name <- gsub("outerroot","roots",coord$Sample.name)

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M2Lathyris$Sample.name))] # is any sample missing?

mzcoord <- merge(M2Lathyris,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="/Users/madeleineernst/Documents/PostDoc/NaturePlants/Science/3DModels_MZmine/Final/Coordinates/ELathyris_Model2_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Hirta Model 1    #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M1Hirta <- ClApp[ClApp$ModelSpecies=="M1hirta",]
M1Hirta$PlantPartID <- as.character(M1Hirta$PlantPartID)

ex <- c("fruitextra","fruit1","fruit2","fruit3","fruit4","fruit5","fruit6")
wh <- which(M1Hirta$PlantPartID %in% ex)  
M1Hirta$PlantPartID[wh] <- c("fruit1","fruit2","fruit3","fruit4","fruit5","fruit6","fruit7")[match(M1Hirta$PlantPartID[wh],ex)]
M1Hirta$PlantPartID <- as.factor(M1Hirta$PlantPartID)

colnames(M1Hirta)[colnames(M1Hirta)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EHirta_coordinates.csv", sep=";", dec=".",header = TRUE) 

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M1Hirta$Sample.name))] # is any sample missing?

mzcoord <- merge(M1Hirta,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="EHirta_Model1_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Hirta Model 2    #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M2Hirta <- ClApp[ClApp$ModelSpecies=="M2hirta",]
colnames(M2Hirta)[colnames(M2Hirta)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EHirta_coordinates.csv", sep=";", dec=".",header = TRUE) 

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M2Hirta$Sample.name))] # is any sample missing?

mzcoord <- merge(M2Hirta,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]
mzcoord2 <- mzcoord2[-which(mzcoord2$Sample.name%in%notin),]

write.table(mzcoord2,file="EHirta_Model2_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Horrida Model 1  #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M1Horrida <- ClApp[ClApp$ModelSpecies=="M1horrida",]
colnames(M1Horrida)[colnames(M1Horrida)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EHorrida_Model1_coordinates.csv", sep=";", dec=".",header = TRUE) 
coord$Sample.name <- gsub("innerroot","roots6",coord$Sample.name)
coord$Sample.name <- gsub("outerroot","roots6",coord$Sample.name)
coord$Sample.name <- gsub("fruits","fruit1",coord$Sample.name)

# 135_Z1bfruits_horridaM1_P1.C.6_01_36459.mzXML is pool2, was labeled wrongly
M1Horrida$Sample.name[M1Horrida$FileName=="135_Z1bfruits_horridaM1_P1.C.6_01_36459.mzXML"] <- "pool2"

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M1Horrida$Sample.name))] # is any sample missing?

mzcoord <- merge(M1Horrida,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="EHorrida_Model1_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Horrida Model 2  #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M2Horrida <- ClApp[ClApp$ModelSpecies=="M2horrida",]

M2Horrida$PlantPartID <- gsub("roots7","innerroot",M2Horrida$PlantPartID)
M2Horrida$PlantPartID <- gsub("outerroot8","outerroot",M2Horrida$PlantPartID)
M2Horrida$PlantPartID <- gsub("flowers1","flowers",M2Horrida$PlantPartID)
colnames(M2Horrida)[colnames(M2Horrida)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EHorrida_Model2_coordinates.csv", sep=";", dec=".",header = TRUE) 

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M2Horrida$Sample.name))] # is any sample missing?

mzcoord <- merge(M2Horrida,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="EHorrida_Model2_features_MZmine.csv",row.names=F,sep=",",quote=T)

####################
# Horrida Model 3  #
####################

ClApp <- read.csv("FeatureTable_MetaData_3DMolecularCartography.csv",sep=",",header = T)
M3Horrida <- ClApp[ClApp$ModelSpecies=="M3horrida",]
M3Horrida$PlantPartID <- gsub("innerroot7","innerroot",M3Horrida$PlantPartID)
M3Horrida$PlantPartID <- gsub("outerroot8","outerroot",M3Horrida$PlantPartID)
M3Horrida$PlantPartID <- gsub("stems1","pool1",M3Horrida$PlantPartID)
colnames(M3Horrida)[colnames(M3Horrida)=="PlantPartID"] <- "Sample.name"

# read coordinates
coord <-read.csv("EHorrida_Model3_coordinates.csv", sep=";", dec=".",header = TRUE) 

notin <- unique(coord$Sample.name)[-which(unique(coord$Sample.name) %in% unique(M3Horrida$Sample.name))] # is any sample missing?

mzcoord <- merge(M3Horrida,coord,by="Sample.name",all=T,sort=F)

# order data frame (Sample.name, x, y, z, r, features)
Sn <- which(colnames(mzcoord) == "Sample.name")
x <- which(colnames(mzcoord) == "x")
y <- which(colnames(mzcoord) == "y")
z <- which(colnames(mzcoord) == "z")
r <- which(colnames(mzcoord) == "r")
features <- which(startsWith(colnames(mzcoord),"X")==T)

mzcoord2 <- mzcoord[,c(Sn,x,y,z,r,features)]

write.table(mzcoord2,file="EHorrida_Model3_features_MZmine.csv",row.names=F,sep=",",quote=T)
