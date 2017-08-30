## Script to format PhormMacroInvert 2016 cyanotoxin and dry wt data

format_PhormMacroInvert2016 <- function(){

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
  inv.wt <-  inv.wt[match(unique(inv.wt$UCSC_ID), inv.wt$UCSC_ID), -c(8:9)]

## Format columns and sample IDs
  inv.wt$date.collected <- as.Date(inv.wt$date.collected, format= "%d-%b-%y")
  inv.wt$site <- revalue(inv.wt$site, c("SFCed"="sfced", "ChimTree"="chimtree", "Pville_DS"="pvilleDS", "USGS_DS"="usgsDS"))
  inv.wt$treatment <- revalue(inv.wt$treatment, c("Cob"="cob", "Pcob"="pcob", "Ph"="phorm"))
  names(inv.wt)[c(6, 8)] <- c("rep", "weight_g")
  inv.wt$sampID <- with(inv.wt, paste(site, treatment, rep, sep="_"))


#### Phormidium and Cladophora weights
  substrate.wt.orig <- as.data.frame(toxin.wt.data$substrate.wt)
  phorm.clad.wt <- subset(substrate.wt.orig, treatment== "phorm_mat" | treatment== "clado")
  phorm.clad.wt$date.collected <- as.Date(phorm.clad.wt$date.collected, format= "%d-%b-%y")
  phorm.clad.wt <- droplevels(phorm.clad.wt)

## Format sample IDs
  phorm.clad.wt$site <- gsub("_.*", "", phorm.clad.wt$ID)
  phorm.clad.wt$site <- as.factor(gsub("pville", "pvilleDS", phorm.clad.wt$site))
  phorm.clad.wt$rep <- as.factor(gsub("[^0-9]", "", phorm.clad.wt$ID))
  phorm.clad.wt$treatment <- revalue(phorm.clad.wt$treatment, c("phorm_mat"="phorm"))


#### Format toxin data ####
  tox.orig <- as.data.frame(toxin.wt.data$all.tox)

## Invertebrate toxin data
  inv.tox <- tox.orig[grep("INVERT", tox.orig$UCSC_ID), ]
  names(inv.tox)[1] <- "UCSC_ID_orig"

## Re-format Inv.Toxin UCSC_ID column
 inv.tox$UCSC_ID <- as.character(inv.tox$UCSC_ID_orig)
 inv.tox$UCSC_ID <- gsub("2016.", "", inv.tox$UCSC_ID)
 inv.tox$UCSC_ID <- as.factor(gsub("[^0-9]", "", inv.tox$UCSC_ID))

## Phormidium, Clado, H2O, and Periphyton toxin data
 tox.non.invert <- tox.orig[-grep("INVERT", tox.orig$UCSC_ID), ]


#### Merge toxin and weight data ####
 phorm.clad.df <- merge(phorm.clad.wt, tox.non.invert, by= "UCSC_ID", all.x= T)
 inv.df <- merge(inv.wt, inv.tox, by= "UCSC_ID")

## Combine into a list
 output.list <- list(inv.df=inv.df, phorm.clad.df=phorm.clad.df)

## Standardize toxins by grams of dry weight
 tox.dw <- lapply(output.list, function(x) {
   sapply(x[, 11:16], function(y) y/x$weight_g)
 })
 output.list.dw <- Map(cbind, output.list, tox.dw) # combine lists

## Rename column headers
 output.list.dw <- lapply(output.list.dw, function(x) {
   setNames(x, nm= c(names(x)[1:16], paste(names(x)[17:22], "_dw", sep="" )))
 })

## Sum up MCY congeners
 totMCY.list <- lapply(output.list.dw, function(x) {
   rowSums(x[ ,19:22])
 })
 output.list.dw <- Map(cbind, output.list.dw, "totMCY_dw"= totMCY.list)

 return(output.list.dw)

}