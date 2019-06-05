#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(dplyr)
library(tidyr)

lib <- read.table(args[1], header=TRUE, stringsAsFactors=FALSE, sep='\t')
colnames(lib) <- c('guideid', 'gene', 'guide', 'barcode')

count <- read.table(args[2], header=FALSE, stringsAsFactors=FALSE, sep='\t')

colnames(count) <- c('guide', 'b1', 'b2', 'count')

count$b2in <- paste(
    count$guide, count$b2, sep='.'
) %in% paste(
           lib$guide, lib$barcode, sep='.'
       )

count$barcode <- ifelse(
    count$b2in, count$b2, count$b1
)

cumcount <- count %>% group_by(
                          guide, barcode
                      ) %>% summarize(count=sum(count))

lib <- left_join(
    lib,
    cumcount[c('guide', 'barcode', 'count')],
    by=c('guide', 'barcode')
)

lib <- replace_na(lib, list('count'=0))

write.table(
    lib, args[3], sep='\t', row.names=FALSE, quote=FALSE
)
