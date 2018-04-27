##################################################
#                                                #
#    R script used for creating Figure 4A        #
#                                                #
##################################################

# Occurences of Euphorbia species investigated chemically and Euphorbia-feeding Hyles retrieved
# from GBIF and manually restricted to native areas

## load libraries
library(rgbif)
library(ggplot2)
library(countrycode)

############# retrieve GBIF data
#https://rdrr.io/cran/rgbif/man/gbifmap.html

source("GBIFMAP_CUSTOMIZED.R")

athymalus <- read.csv("Athymalus_0003633-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.ug6cnm
chamaescyce <- read.csv("Chamaesyce_0003634-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.9ntc6n
esula <- read.csv("Esula_0003635-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.avy89i
euphorbia <- read.csv("Euphorbia_0003638-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.llmwkb

# AS = Asia
# EU = Europe
# NA = North America
# OC = Oceania
# AF = Africa
# SA = South America
# AN = Antarctica

dict <- read.csv("country_continent.csv",stringsAsFactors = F,na.strings = "NNN")

###############
## Athymalus  #
###############
# horrida has no data
athymalus_nat <- athymalus
ath_conts <- countrycode(unique(athymalus_nat$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
athymalus_nat$continent <- ath_conts[match(athymalus_nat$countrycode,unique(athymalus_nat$countrycode))]
athymalus_nat <- athymalus_nat[athymalus_nat$continent=="AF"|athymalus_nat$continent=="AS",]

colnames(athymalus_nat)[colnames(athymalus_nat)=="decimallatitude"] <- "decimalLatitude"
colnames(athymalus_nat)[colnames(athymalus_nat)=="decimallongitude"] <- "decimalLongitude"
colnames(athymalus_nat)[colnames(athymalus_nat)=="species"] <- "name"

gbifmap_custom(athymalus_nat)

################
## Chamaesyce  #
################

chamaesyce_nat <- chamaescyce
cham_conts <- countrycode(unique(chamaesyce_nat$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
chamaesyce_nat$continent <- cham_conts[match(chamaesyce_nat$countrycode,unique(chamaesyce_nat$countrycode))]
chamaesyce_nat <- chamaesyce_nat[-which(chamaesyce_nat$species=="Euphorbia thymifolia" & !chamaesyce_nat$continent %in% c("SA","AS","NA")),]
chamaesyce_nat <- chamaesyce_nat[-which(chamaesyce_nat$species != "Euphorbia platyclada" & chamaesyce_nat$continent == "AF"),]
chamaesyce_nat <- chamaesyce_nat[-which(chamaesyce_nat$species != "Euphorbia platyclada" & chamaesyce_nat$species != "Euphorbia thymifolia" & !chamaesyce_nat$continent %in% c("SA","NA")),]
chamaesyce_nat <- chamaesyce_nat[-which(chamaesyce_nat$species == "Euphorbia platyclada" & chamaesyce_nat$continent == "NA"),]

colnames(chamaesyce_nat)[colnames(chamaesyce_nat)=="decimallatitude"] <- "decimalLatitude"
colnames(chamaesyce_nat)[colnames(chamaesyce_nat)=="decimallongitude"] <- "decimalLongitude"
colnames(chamaesyce_nat)[colnames(chamaesyce_nat)=="species"] <- "name"

gbifmap_custom(chamaesyce_nat)

###########
## Esula  #
###########
# acanthothamnos has no data
esula_nat <- esula
esula_conts <- countrycode(unique(esula_nat$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
esula_nat$continent <- esula_conts[match(esula_nat$countrycode,unique(esula_nat$countrycode))]
esula_nat <- esula_nat[esula_nat$continent=="AS"|esula_nat$continent=="EU",]
esula_nat <- esula_nat[-which(esula_nat$species == "Euphorbia amygdaloides" & esula_nat$continent=="AS"),]
esula_nat <- esula_nat[-which(esula_nat$species == "Euphorbia characias" & esula_nat$continent != "EU"),]
esula_nat <- esula_nat[-which(esula_nat$species == "Euphorbia cyparissias" & esula_nat$continent == "AS"),]
esula_nat <- esula_nat[-which(esula_nat$species == "Euphorbia lathyris" & esula_nat$continent == "EU"),]

colnames(esula_nat)[colnames(esula_nat)=="decimallatitude"] <- "decimalLatitude"
colnames(esula_nat)[colnames(esula_nat)=="decimallongitude"] <- "decimalLongitude"
colnames(esula_nat)[colnames(esula_nat)=="species"] <- "name"

gbifmap_custom(esula_nat)

###############
## Euphorbia  #
###############
# abdelkuri has no data
# lactea has no data in native range
# stenoclada has one record in Pacific Ocean, can be removed manually
euphorbia_nat <- euphorbia
euphorbia_conts <- countrycode(unique(euphorbia_nat$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
euphorbia_nat$continent <- euphorbia_conts[match(euphorbia_nat$countrycode,unique(euphorbia_nat$countrycode))]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia alluaudii" & euphorbia_nat$continent=="EU"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia cylindrifolia" & euphorbia_nat$continent=="EU"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia grandicornis" & euphorbia_nat$continent=="NA"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia ingens" & euphorbia_nat$continent %in% c("SA","EU","NA","OC")),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia lactea"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia milii" & euphorbia_nat$countrycode!="MG"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia neriifolia" & euphorbia_nat$continent %in% c("SA","NA","EU")),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia sipolisii" & euphorbia_nat$continent=="EU"),]
euphorbia_nat <- euphorbia_nat[-which(euphorbia_nat$species=="Euphorbia stenoclada" & euphorbia_nat$continent %in% c("EU","OC")),]

colnames(euphorbia_nat)[colnames(euphorbia_nat)=="decimallatitude"] <- "decimalLatitude"
colnames(euphorbia_nat)[colnames(euphorbia_nat)=="decimallongitude"] <- "decimalLongitude"
colnames(euphorbia_nat)[colnames(euphorbia_nat)=="species"] <- "name"

gbifmap_custom(euphorbia_nat)

######################
## Hyles euphorbiae  #
######################

hyles <- read.csv("Hyles_0003640-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.cgwykr
hyles_conts <- countrycode(unique(hyles$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
hyles$continent <- hyles_conts[match(hyles$countrycode,unique(hyles$countrycode))]
hyles <- hyles[hyles$year<1990,]

colnames(hyles)[colnames(hyles)=="decimallatitude"] <- "decimalLatitude"
colnames(hyles)[colnames(hyles)=="decimallongitude"] <- "decimalLongitude"
colnames(hyles)[colnames(hyles)=="species"] <- "name"

##################################
## Other Euphorbia-feeding Hyles #
##################################
#[Hundsdoerfer et al., 2009]

efhyles <- read.csv("EFeedingHyles_0012430-171219132708484.csv", sep="\t",stringsAsFactors = F,quote="") # downloaded from https://doi.org/10.15468/dl.ahwr3w
efhyles_conts <- countrycode(unique(efhyles$countrycode),custom_dict= dict, origin = "iso.3166.country",destination="continent.code")
efhyles$continent <- efhyles_conts[match(efhyles$countrycode,unique(efhyles$countrycode))]

colnames(efhyles)[colnames(efhyles)=="decimallatitude"] <- "decimalLatitude"
colnames(efhyles)[colnames(efhyles)=="decimallongitude"] <- "decimalLongitude"
colnames(efhyles)[colnames(efhyles)=="species"] <- "name"

############ plot all

all <- rbind(athymalus_nat, chamaesyce_nat,esula_nat,euphorbia_nat)
all <- rbind(all, hyles,efhyles)

subgenera <- read.csv("Euphorbia_Subgenera.csv",stringsAsFactors = F,header = F)
all$subgenera <- subgenera$V2[match(all$name, subgenera$V1)]
all <- all[!is.na(all$subgenera),] # will remove E. polygona (can be disregarded)

colnames(all)[colnames(all)=="name"] <- "species"
colnames(all)[colnames(all)=="subgenera"] <- "name"

gbifmap_custom(all)

