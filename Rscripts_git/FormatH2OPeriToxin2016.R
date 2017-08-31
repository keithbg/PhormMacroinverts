## Format H2O, Periphyton, and Rinse cyanotoxin data collect as part of
## Phormidium Macroinvertebrates dataset in 2016


## Read in import function, which outputs a list of 3 data frames
  source("/Users/KeithBG/Documents/UC Berkeley/Manuscripts/PhormMacroinverts/Rscripts_git/ImportInvertToxin2016.R")
  toxin.wt.data <- import_PhormMacroInvert_data()

#### H2O, periphyton, and rinse (hpr) weights and ID formating ####
  substrate.wt.orig <- as.data.frame(toxin.wt.data$substrate.wt)
  hpr.wt <- subset(substrate.wt.orig, dataset!= "phorm_inverts")
  hpr.wt$date.collected <- as.Date(hpr.wt$date.collected, format= "%d-%b-%y")
  hpr.wt <- droplevels(hpr.wt)

## Add area column for area of periphyton scraped (cm^2)
  # 62.2 cm^2 =  area under a 1L Nalgene cap
  hpr.wt$area_cm2 <- as.numeric(NA)
  hpr.wt[hpr.wt$dataset== "periphyton", "area_cm2"] <- 62.2
  hpr.wt[hpr.wt$date.collected == "2016-06-05","area_cm2"] <- 62.2*3

## Add 1L of volume filtered for H2O samples
  h2o.vol <- read.table("/Users/KeithBG/Documents/UC Berkeley/2016 Summer Research/PhormidiumMacroinverts2016/h2oFilterVol.tab", sep= "\t", header= TRUE)
  hpr.wt[hpr.wt$dataset== "H2O", "vol_L"] <- h2o.vol$vol_L

## Format site and rep column
  hpr.wt$site <- gsub("_.*", "", hpr.wt$treatment)
  hpr.wt$site <- as.factor(gsub("pville", "pvilleDS", hpr.wt$site))
  hpr.wt$rep <- as.factor(gsub("[^0-9]", "", hpr.wt$treatment))


#### Format cyanotoxin data #####



#### Merge weights/volumes/area and cyanotoxin data ####




