#################################################
#                                               #
#    R script used for creating Figure S10      #
#                                               #
#################################################

# Differential expression (TIC normalized MS1 intensities) of molecules annotated as Euphorbia diterpenoids
# in different plant parts of one representative species per Euphorbia subgeneric clade

## load library
library(pheatmap)

## read mass spectral molecular network data
net <-read.csv("Cytoscape_SummaryTable_3DMolecularCartography.csv", sep=",", dec=".",header = TRUE) # Table exported from Cytoscape, with all network attributes loaded

IDs <- paste(round(net$row.m.z,4),round(net$row.retention.time,4),paste("ID",net$shared.name,sep=".."),sep="_")
net <- cbind(IDs, net)

# filter out only feature IDs annotated as diterpenoids
c1 <- which(colnames(net) %in% c("LibraryID","IDs","componentindex","Diterpenoid_Subclass_Interpretation"))
c2 <- grep("Peak.area",colnames(net))

net <- net[,c(c1,c2)]

dit <- c("DSF 313, 295, 285 Euphorbia diterpenoid",
         "DSF 313, 295, Euphorbia diterpenoid",
         "DSF 317, 299, 291, 281 Euphorbia diterpenoid",
         "DSF 329, 311, 293, 283 - Euphorbia diterpenoid",
         "DSF 337, 313, 295, 277, 267 Euphorbia diterpenoid",
         "Euphorbia factor L14 - Lathyrane diterpenoid",
         "Euphorbia factor L8 - Lathyrane diterpenoid - likely",
         "Ingenol O-Decadienoyl, O-Ac or isomer putative",
         "Miliiamin C derivative",
         "Miliiamin C new derivative",
         "Miliiamine C 20-Ac, (or positional isomer)",
         "Miliiamine C 20-Ac, N-de-Me - Ingenol diterpenoid",
         "Milliamine C; 5-Anthraniloyl type moiety, 20-Ac - Ingenol diterpenoid",
         "m/z 295 of diterpene esters (collected from 12-deoxyphorbol-13-acetate)",
         "m/z 313 of diterpene esters (collected from 12-deoxyphorbol-13-acetate)")
          
net <- net[net$LibraryID %in% dit | net$Diterpenoid_Subclass_Interpretation =="Euphorbia diterpenoids",]

ClApp <- net

# transpose and order to plot heatmap
ClApp <- t(ClApp)
colnames(ClApp) <- ClApp[1,]
ClApp <- ClApp[-c(1:4),]
class(ClApp) <- "numeric"

# are there any 0 columns?
which(colSums(ClApp)==0)

# subsitute wrong columnname
rownames(ClApp)[rownames(ClApp)=="X135_Z1bfruits_horridaM1_P1.C.6_01_36459.mzML.Peak.area"]   <- "X135_Z2apool_horridaM1_P1.C.6_01_36459.mzML.Peak.area"
l <- strsplit(rownames(ClApp),"_")
newr <- paste(sapply(l, "[[", 2),sapply(l, "[[", 3),sep="_")
rownames(ClApp) <- newr

hirta <- grep("hirta",rownames(ClApp))
horrida <- grep("horrida",rownames(ClApp))
lathyris <- grep("lathyris",rownames(ClApp))
milii <- grep("milii",rownames(ClApp))

# hirta
h1 <- grep("leaves",rownames(ClApp)[hirta])
h2 <- grep("stem",rownames(ClApp)[hirta])
h3 <- grep("fruit",rownames(ClApp)[hirta])
h4 <- grep("innerroot",rownames(ClApp)[hirta])
h5 <- grep("outerroot",rownames(ClApp)[hirta])

# horrida
ho1 <- grep("pool",rownames(ClApp)[horrida])
ho2 <- grep("spines",rownames(ClApp)[horrida])
ho3 <- grep("flowers",rownames(ClApp)[horrida])
ho4 <- grep("fruits",rownames(ClApp)[horrida])
ho5 <- grep("root",rownames(ClApp)[horrida])

# lathyris
l1 <- grep("leaves",rownames(ClApp)[lathyris])
l2 <- grep("stems",rownames(ClApp)[lathyris])
l3 <- grep("root",rownames(ClApp)[lathyris])

# milii
m1 <- grep("leaves",rownames(ClApp)[milii])
m2 <- grep("stems",rownames(ClApp)[milii])
m3 <- grep("thorns",rownames(ClApp)[milii])
m4 <- grep("flowers",rownames(ClApp)[milii])
m5 <- grep("roots",rownames(ClApp)[milii])

# sort rows per species
df <- ClApp[c(hirta[h1],hirta[h2],hirta[h3],hirta[h4],hirta[h5],
              horrida[ho1],horrida[ho2],horrida[ho3],horrida[ho4],horrida[ho5],
              lathyris[l1],lathyris[l2],lathyris[l3],
              milii[m1],milii[m2],milii[m3],milii[m4],milii[m5]),]

makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

cutoff.distance <- 0.001
cols <- makeColorRampPalette(c("#2C7BB6", "yellow",   
                               "yellow", "red"), 
                             cutoff.distance / max(df),ncol(df))

pdf("Heatmap_3DMolecularCartography.pdf",width=8,height=5, bg="transparent")
pheatmap(df,  clustering_method="average", clustering_distance_cols="maximum", 
         show_colnames=TRUE,show_rownames=TRUE, cluster_rows=FALSE, cluster_cols=TRUE, 
         fontsize_row=2,fontsize_col=1,color=cols)
dev.off()

