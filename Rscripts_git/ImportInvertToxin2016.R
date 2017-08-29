## Function to import and format cyanotoxin concentrations and column names
## from Phorm invert samples collected 2016.
## Function returns a list with 3 elements: 1) all toxin data, 2) invert dry weights
## and 3) dry weights and volumes of all substrate scrapings


## Cyanotoxin data from Raphe Kudela's LCMS units ng toxin / mL
## Samples processed in UCSC on 3/4-Nov-2016

## Function pathway:
## /Users/KeithBG/Documents/UC Berkeley/Manuscripts/PhormMacroinverts/Rscripts_git/")

import_PhormMacroInvert_data <- function(){

## Location of 2016 data
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
## Multiply MCY concentrations by 4 to account for subsampling
  # 6 mL total; 3 mL subsample for SPE cleaning;
  # 2mL produced after SPE cleaning; 1 mL subsample put into LC-MS vial
  # (6/3) * (2/1) = 4
  cob.tox[, which(names(cob.tox) %in% c("LR", "RR", "YR", "LA"))] <- 4*cob.tox[, which(names(cob.tox) %in% c("LR", "RR", "YR", "LA"))]

## Multiply ATX concentration by 6 to account for subsampling
  # 6mL total with 1mL subsample into LC-MS vial: 6/1 = 6
  cob.tox[, which(names(cob.tox) %in% c("ATX"))] <- 6*cob.tox[, which(names(cob.tox) %in% c("ATX"))]

## Inv samples were extracted in 2ml 50% MeOH and 1 mL subsamples for both ATX and MCY
  # 2/1= 2
  inv.tox[, which(names(inv.tox) %in% c("ATX", "LR", "RR", "YR", "LA"))] <- 2*inv.tox[, which(names(inv.tox) %in% c("ATX", "LR", "RR", "YR", "LA"))]

## Combine inv and cob toxins
  all.tox <- rbind(inv.tox, cob.tox)
  names(all.tox)[1] <- "UCSC_ID"


#### Dry Weight Data ####

#### Macroinvertebrate dry weights
  filename <- paste(pathway_2016, "MacroinvertWeights2016.tab", sep= "/")
  inv.wt <- read.table(filename, header= TRUE, sep= "\t")

## Remove incomplete rows
  inv.wt <- subset(inv.wt, UCSC_ID != "<DL" & UCSC_ID != "")
  inv.wt <- inv.wt[grep("[0-9]", inv.wt$UCSC_ID), ] # only rows with numeric ID's

## Format column names
  names(inv.wt) <- c("date.collected", "site", "treatment", "order", "family", "cob.num", "UCSC_ID", "weight_g", "notes")

#### Phormidium and cobble scraping dry weights
  filename <- paste(pathway_2016, "CyanotoxinParticulateSampleIDs.tab", sep= "/")
  substrate.wt <- read.table(filename, header= TRUE, sep= "\t")

## Combine all data frames into a list
  toxin.wt.list <- list(all.tox= all.tox, inv.wt= inv.wt, substrate.wt= substrate.wt)
  return(toxin.wt.list)

}