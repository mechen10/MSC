#!/usr/bin/env Rscript
# Reorders files for primer


library(optparse)

option_list = list(
  make_option(c("-m", "--mappingfile"), type="character", default=NULL, 
              help="mappingfile fp ", metavar="character"),
  make_option(c("-d", "--dm"), type="character", 
              help="distance matrix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="OUTPUT", 
              help="output name for 2 files", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$mappingfile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (mappingfile).n", call.=FALSE)
}

if (is.null(opt$dm)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (dm).n", call.=FALSE)
}

metadata <- read.delim(opt$mappingfile
                       , header = TRUE
                       , row.names = 1
                       , na.strings = c("","NA","na")
                       , stringsAsFactors = FALSE
                       , strip.white = TRUE)

dm <- read.delim(opt$dm
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE
                 , strip.white = TRUE)
# For some reason, they change headers to '.' so I'm changing them back
colnames(dm) <- gsub(".", "-", colnames(dm), fixed =TRUE)

# Find duplicate values and get rid of missing ones
combined <- c(row.names(metadata), row.names(dm))
combined.dup <- combined[duplicated(combined)]

positionsMetadata <- lapply(combined.dup, function(x) grep(x, rownames(metadata)))
positionsMetadata <- unlist(positionsMetadata)
metadataEdited <- metadata[positionsMetadata,]

positionsDMrow <- lapply(combined.dup, function(x) grep(x, rownames(dm)))
positionsDMcol <- lapply(combined.dup, function(x) grep(x, colnames(dm)))
positionsDMrow <- unlist(positionsDMrow)
positionsDMcol <- unlist(positionsDMcol)

dmEdited <- dm[positionsDMrow,positionsDMcol]

# Now, reorder them alphabetically

metadataEdited <- metadataEdited[order(rownames(metadataEdited)),]

dmEdited <- dmEdited[order(rownames(dmEdited)),]
dmEdited <- dmEdited[,order(colnames(dmEdited))]
names <- rownames(dmEdited)
dmEditedPrimer <- cbind(names,dmEdited)


# Finally, insert 2 columns for mapping file

# Make Factor
Factor <- rep(1,length(rownames(metadataEdited)))
BLANK <- rep("",length(rownames(metadataEdited)))
dfInsert <- cbind(Factor,BLANK)
dfInsert <- as.data.frame(dfInsert)
rownames(dfInsert) <- rownames(metadataEdited)

metadataEditedPrimer <- merge(dfInsert, metadataEdited, by = 0)
colnames(metadataEditedPrimer)[1] <- "#SampleID"

colnames(metadataEditedPrimer)[3] <- ""
colnames(dmEditedPrimer)[1] <- ""


write.table(metadataEditedPrimer, file = paste0(opt$out,"-metadata.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(dmEditedPrimer, file = paste0(opt$out,"-dm.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

