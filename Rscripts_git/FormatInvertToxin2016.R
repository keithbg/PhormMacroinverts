## Script to format PhormMacroInvert 2016 cyanotoxin and dry wt data

library(plyr)

## Read in import function, which outputs a list of 3 data frames
  source("/Users/KeithBG/Documents/UC Berkeley/Manuscripts/PhormMacroinverts/Rscripts_git/ImportInvertToxin2016.R")
  toxin.wt.data <- import_PhormMacroInvert_data()

#### Format wt data ####

#### Invertebrate weights
  inv.wt.orig <- as.data.frame(toxin.wt.data$inv.wt)
  inv.wt.orig <- droplevels(inv.wt.orig)
  inv.wt.orig$weight_g <- as.numeric(as.character(inv.wt.orig$weight_g))

## Combine samples that went into same UCSC vial
  inv.wt <- ddply(inv.wt.orig, .(UCSC_ID), mutate,
                samp.wt = sum(weight_g))
  inv.wt <-  test[match(unique(test$UCSC_ID), test$UCSC_ID), -c(8:9)]

## Format columns
  inv.wt$sampID <- with(inv.wt, paste(site, treatment, cob.num, sep="_"))
  inv.wt$date.collected <- as.Date(inv.wt$date.collected, format= "%d-%b-%y")

#### Phormidium and Cladophora weights
  substrate.wt.orig <- as.data.frame(toxin.wt.data$substrate.wt)
  phorm.clad.wt <- subset(substrate.wt.orig, treatment== "phorm_mat" | treatment== "clado")

#### Format toxin data ####
  tox.orig <- as.data.frame(toxin.wt.data$all.tox)

#### Invertebrate toxin data
  tox.inv <- tox.orig[grep("INVERT", tox.orig$UCSC_ID), ]
  names(tox.inv)[1] <- "UCSC_ID_orig"

## Re-format Inv.Toxin UCSC_ID column
 tox.inv$UCSC_ID <- as.character(tox.inv$UCSC_ID_orig)
 tox.inv$UCSC_ID <- gsub("2016.", "", tox.inv$UCSC_ID)
 tox.inv$UCSC_ID <- as.factor(gsub("[^0-9]", "", tox.inv$UCSC_ID))

#### Phormidium, Clado, H2O, and Periphyton toxin data
 tox.other <- tox.orig[-grep("INVERT", tox.orig$UCSC_ID), ]


#### Merge toxin and weight data ####
 phorm.clad.df <- merge(phorm.clad.wt, tox.other, by= "UCSC_ID", all.x= T)


