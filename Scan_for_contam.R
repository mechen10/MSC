#!/bin/bash
### Scan for contaminants ###
###################### DESCRPTION OF SCRIPT #######################
# This script scans the whole OTU table using the control sequences and figures out if there are contaminants that we should remove.
# Typical characteristics of a contaminant are:
  # It is present in the blank or control sequences (duh!)
  # It is relatively abundant in low-read samples and relatively spares in high-read samples

# Thus, this script takes the raw OTU table and plots suspected contaminants.
# It does NOT remove OTUs from the table; it simply helps you identify potential contaminants so you can further investigate them.

library(optparse)

option_list <- list(
  make_option(c("-i", "--otu_table"), default = "EMPTY", help="OTU table (raw) in text or biom format."),
  make_option(c("-m", "--mapping_file"), default = "EMPTY", help="Mapping file or metadata."),
  make_option(c("-b", "--blanks"), help="The column and treatment for blanks. Eg.'Cateogry:treatment1,treatment2'"),
  make_option(c("-s", "--split_table"), help="If you want to split table into multiple tables, indicate column and table to include for each table from mapping file. Eg. 'Column:table1,table2'"),
  make_option(c("--filter_low_abund_per_sample"), default=0, help="Prior to scanning, change all OTUs in each sample below this threshold to zero. Default: 0"),
  make_option(c("--filter_low_abund_overall"), default=0, help="Prior to scanning, remove all OTUs whose maximum read count in any sample is less than this threshold. (Note that this is done BEFORE filtering by sample) Default: 0"),
  make_option(c("-o", "--output"), default="control_analysis", help="Folder name for output files"),
  make_option(c("--dashes"), default="FALSE", help="If your sample names have dashes, mark this as TRUE"),
  make_option(c("--numberstart", default="FALSE", help="If your sampe names start with numbers, mark this as TRUE"))
  )
opt = parse_args(OptionParser(option_list=option_list))
otu_table = opt$otu_table
mf_fp = opt$mapping_file
blanks = opt$blanks
splittable = opt$split_table
thresh_sample = opt$filter_low_abund_per_sample
thresh_all = opt$filter_low_abund_overall
output = opt$output
dashes = opt$dashes
numberstart = opt$numberstart

# Exit if something is missing
if (otu_table == "EMPTY") {
  print("Please provide OTU Table")
  quit()
}
if (mf_fp == "EMPTY") {
  print("Please provide MF")
  quit()
}

############# FOR TESTING###################
# setwd("/Users/mckenzielab/Documents/Melissa/Project_AM/control_analysis/")
# otu_table <- "/Users/mckenzielab/Documents/Melissa/Project_AM/OTU_Table_wtaxa.biom"
# mf_fp <- "/Users/mckenzielab/Documents/Melissa/Project_AM/Merged_mapping_file.txt"
# blanks <- "Description:CON"
# splittable <- "Project:H,P"
# thresh_sample <- "5"
# thresh_all <- "10"
# output <- "TESTING"
# dashes <- TRUE
# numberstart <- TRUE



################# READ IN DATA ###############
print("Reading in data...")

# Make output folder
dir.create(output)
setwd(output)
system(paste0("biom convert -i ",otu_table," --to-tsv --header-key taxonomy -o OTU_Table_text_forcontamscan.txt"))

otu <- read.delim("OTU_Table_text_forcontamscan.txt", header=TRUE, row.names=1, skip =1, stringsAsFactors = FALSE)
mf <- read.delim(paste0(mf_fp), header=TRUE, row.names=1, stringsAsFactors = FALSE, na.strings = "")

# Get blank names
blanks_split <- unlist(strsplit(blanks, split = c(":")))
blanks_split <- c(blanks_split[1], strsplit(blanks_split[2], ","))

# Get items to split table by
table_split <- unlist(strsplit(splittable, split = c(":")))
table_split <- c(table_split[1], strsplit(table_split[2], ","))


###################### SETTING UP TABLES #########################
print("Setting up tables...")
# Get rid of taxonomy
taxa <- data.frame(rownames(otu), otu[,ncol(otu)])
otu <- otu[,-ncol(otu)]
# Fix header names if applicable
if (numberstart == "TRUE") {
  colnames(otu) <- gsub("^X","",colnames(otu))
}
if (dashes == "TRUE") {
  rownames(mf) <- gsub("-",".", rownames(mf), fixed=TRUE)
}

# Check to make sure mf and otu are the same; first, make sure they are shared

tofilt <- which(rownames(mf) %in% colnames(otu)) 
mf.filt <- mf[tofilt,]

# Not, split up into different studies
otu.bystudy <- list()
mf.bystudy <- list()
for ( study in table_split[[2]] ) {
  mf.bystudy[[study]] <- mf.filt[which(mf.filt$Project == study),]
  otu.bystudy[[study]] <- otu[,match(rownames(mf.bystudy[[study]]), colnames(otu))]
}

###### (1) CHANGE ALL COUNTS WITH LESS THAN 5 PER SAMP TO 0 ######

for ( study in names(otu.bystudy) ) {
  otu.temp <- otu.bystudy[[study]]
  for ( r in 1:nrow(otu.temp) ) {
    for (c in 1:ncol(otu.temp) ) {
      if (otu.temp[r,c] <5) {
        otu.temp[r,c] <- 0
      }
    }
  }
  otu.bystudy[[study]] <- otu.temp
}

###### (2) REMOVE ALL OTUS WITH LESS THAN 5 IN WHOLE OTU TABLE ####
# remove all zeros from otu table

for ( study in names(otu.bystudy) ) {
  otu.temp <- otu.bystudy[[study]]
  todel <- c()
  for ( r in 1:nrow(otu.temp)) {
    if (max(otu.temp[r,]) < 5) {
      todel <- c(todel, r)
    }
  }
  otu.bystudy[[study]] <- otu.temp[-todel,]
}


###### (3) FILTER MF AND OTU TO MAKE SURE ALL SAMPLES INCLUDED ######
# Check to make sure everything is in everything
for ( study in names(otu.bystudy) ) {
  # all otu in mf?
  print(any(!(colnames(otu.bystudy[[study]]) %in% rownames(mf.bystudy[[study]]))))
  # all mf in otu?
  print(any(!(rownames(mf.bystudy[[study]]) %in% colnames(otu.bystudy[[study]]))))
}

for ( study in names(otu.bystudy) ) {
  mf.temp <- mf.bystudy[[study]]
  otu.temp <- otu.bystudy[[study]]
  # Get sum of all sample reads, and also filter based on which are actually present in the otu table
  contam.temp <- otu.temp[,mf.temp$Description=="CON"]
  names.contam.temp <- rownames(contam.temp[rowSums(contam.temp) >0,])
  otu.temp.filt.nocon <- otu.temp[,!c(mf.temp$Description=="CON")]
  RPS.temp <- colSums(otu.temp.filt.nocon)
  otu.temp.filt <- otu.temp.filt.nocon[names.contam.temp,]
  # make list of lm of reads in sample vs size of sample
  all.lm.temp <- list()
  for ( r in 1:nrow(otu.temp.filt)) {
    temp <- lm(as.numeric(otu.temp.filt[r,])/RPS.temp ~ as.vector(RPS.temp))
    m <- coef(temp)[2]
    p <- summary(temp)$coefficients[2,4]
    all.lm.temp[[paste0(rownames(otu.temp.filt)[r])]] <- c(m,p)
  }
  # if the slope is negative AND it's significant, then same as tent contam
  tent.contam.temp <- c()
  not.contam.temp <- c()
  for ( o in names(all.lm.temp) ) {
    temp <- all.lm.temp[[paste0(o)]]
    if ( (temp[2] < 0.05) & (temp[1] < 0) ) {
      tent.contam.temp <- c(tent.contam.temp, o)
    } else {
      not.contam.temp <- c(not.contam.temp, o)
    }
  }
  
  # now plot contam and get names
  names.contam.temp <- taxa[match(tent.contam.temp,taxa[,1]),2]
  split.names.temp <- strsplit(as.character(names.contam.temp), split = "; __")
  # get abundance in controls
  contam.tab <- contam.temp[tent.contam.temp,]
  contam.counts <- c()
  for ( r in 1:nrow(contam.tab) ) {
    toprint <- ""
    for ( c in 1:ncol(contam.tab) ) {
      toprint <- paste0(toprint," ",names(contam.tab)[c],":",contam.tab[r,c])
    }
    contam.counts[r] <- toprint
  }
  mat.dim <- ceiling(sqrt(length(tent.contam.temp))) # make matrix
  pdf(file = paste0("Possible_contam_",study,".pdf"), width = 10*(mat.dim/3), height = 10*(mat.dim/3))
  par(mfrow=c(mat.dim, mat.dim))
  count <- 1
  for (n in tent.contam.temp) {
    plot(as.numeric(otu.temp.filt.nocon[n,])/RPS.temp ~ RPS.temp, ylab="Rel Abund", main=paste0(split.names.temp[[count]][6], "(",n,")"), xlab="Reads/sample", sub=paste0(contam.counts[count]), cex.sub=0.5)
    count <- count + 1
  }
  dev.off()
  
  # plot non-contam and get names
  names.noncontam.temp <- taxa[match(not.contam.temp,taxa[,1]),2]
  split.names.non.temp <- strsplit(as.character(names.noncontam.temp), split = "; __")
  # get abundance in controls
  ncontam.tab <- contam.temp[not.contam.temp,]
  ncontam.counts <- c()
  for ( r in 1:nrow(ncontam.tab) ) {
    toprint <- ""
    for ( c in 1:ncol(ncontam.tab) ) {
      toprint <- paste0(toprint," ",names(ncontam.tab)[c],":",ncontam.tab[r,c])
    }
    ncontam.counts[r] <- toprint
  }
  
  mat.dim <- ceiling(sqrt(length(not.contam.temp)))# make matrix
  pdf(file = paste0("Other_contam_",study,".pdf"), width = 10*(mat.dim/3), height = 10*(mat.dim/3))
  par(mfrow=c(mat.dim, mat.dim))
  count <- 1
  for (n in not.contam.temp) {
    plot(as.numeric(otu.temp.filt.nocon[n,])/RPS.temp ~ RPS.temp, ylab="Rel Abund", main=paste0(split.names.non.temp[[count]][6], "(",n,")"), xlab="Reads/sample",sub=paste0(ncontam.counts[count]), cex.sub=0.5)
    count <- count + 1
  }
  dev.off()
  
}

system("rm OTU_Table_text_forcontamscan.txt")


