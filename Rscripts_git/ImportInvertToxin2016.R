## Script to import and format cyanotoxin data from invert samples collected 2016
## Cyanotoxin data from Raphe Kudela's LCMS units ng toxin / mL
## Samples processed in UCSC on 3/4-Nov-2016


library(plyr)

setwd("/Users/KeithBG/Documents/UC Berkeley/Manuscripts/PhormMacroinverts/Rscripts_git/")
pathway_2016 <- "/Users/KeithBG/Documents/UC Berkeley/2016 Summer Research/PhormidiumMacroinverts2016"
#### Cyanotoxin Data ####

#### Read in macroinvert toxin data
  filename <- paste(pathway_2016, "CyanotoxinsInverts2016.tab", sep= "/")
  inv.tox <- read.table(filename, header= TRUE, sep= "\t")

#### Read in particulate cyanotoxin data
  filename <- paste(pathway_2016, "CyanotoxinsParticulate2016.tab", sep= "/")
  part.tox <- read.table(filename, header= TRUE, sep= "\t")
  part.tox <- part.tox[complete.cases(part.tox), ]

## Subset only the Phormidium and periphyton rows
  cob.tox <- part.tox[c(1:24, 42:53), ]
  rm(part.tox)

#### Format toxin columns
## Multiply MCY concentrations by 4 to account for subsampl=ing
  # 6 mL total; 3 mL subsample for SPE cleaning;
  # 2mL produced after SPE cleaning; 1 mL subsample put into LC-MS vial
  # (6/3) * (2/1) = 4
  inv.tox[, which(names(inv.tox) %in% c("LR", "RR", "YR", "LA"))] <- 4*inv.tox[, which(names(inv.tox) %in% c("LR", "RR", "YR", "LA"))]
  cob.tox[, which(names(cob.tox) %in% c("LR", "RR", "YR", "LA"))] <- 4*cob.tox[, which(names(cob.tox) %in% c("LR", "RR", "YR", "LA"))]

## Multiply ATX concentration by 6 to account for subsampling
  # 6mL total with 1mL subsample into LC-MS vial: 6/1 = 6
  inv.tox[, which(names(inv.tox) %in% c("ATX"))] <- 6*inv.tox[, which(names(inv.tox) %in% c("ATX"))]
  cob.tox[, which(names(cob.tox) %in% c("ATX"))] <- 6*cob.tox[, which(names(cob.tox) %in% c("ATX"))]


#### Dry Weight Data ####

#### Macroinvertebrate dry weights
  filename <- paste(pathway_2016, "MacroinvertWeights2016.tab", sep= "/")
  inv.wt <- read.table(filename, header= TRUE, sep= "\t")

## Remove incomplete rows
  inv.wt <- subset(inv.wt, UCSC_ID != "<DL" & UCSC_ID != "")
  inv.wt <- inv.wt[grep("[0-9]", inv.wt$UCSC_ID), ] # only rows with numeric ID's

