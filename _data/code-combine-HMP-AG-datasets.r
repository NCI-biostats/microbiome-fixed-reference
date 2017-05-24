rm(list = ls())
require(plyr)
require(dplyr)
require(Hotelling)
require(xtable)
require(ggplot2)
source('../code-microbiome-helper-functions.r')


# generate the HMP and AG file (some edits were done by hand, ie. removing the prefixes in the bacteria names in the respective files)
# AG (AG_L2_030317-mm) bacteria columns: 6+
# HMP (yw-HMP_data_082816-5-sites-level2-with-bacteria-names-mm) bacteria columns: 5+

ag <- read.csv('./AG_L2_030317-mm.csv')
ag$id <- 1:dim(ag)[1]

hmp <- read.csv('./_data/yw-HMP_data_082816-5-sites-level2-with-bacteria-names-mm.csv')


ag.names <- names(ag)
hmp.names <- names(hmp)

n.ag <- dim(ag)[1]
n.hmp <- dim(hmp)[1]

ix <- match(hmp.names, ag.names)

ix.hmp <- !is.na(ix)
ix.ag <- na.omit(ix)

hmp.names[ix.hmp]
ag.names[ix.ag]

d.hmp <- hmp[, ix.hmp]
d.ag <- ag[, ix.ag]


d.hmp.ag <- rbind(cbind(project = 'HMP', d.hmp),
                  cbind(project = 'AG', d.ag))

write.csv(d.hmp.ag, file = './_data/mm-HMP_data_082816-5-sites-level2-AG_L2_030317-with-bacteria-names.csv')


# now include non-overlapping columns
ix.hmp.other <- !ix.hmp

ix.ag.incl.tf <- 1:length(ag.names) %in% ix.ag
ix.ag.other.tf <- !ix.ag.incl.tf
ix.ag.other <- (1:length(ag.names))[ix.ag.other.tf]

d11 <- hmp[, ix.hmp.other]
d21 <- matrix(0, ncol = sum(ix.hmp.other), nrow = n.ag)
colnames(d21) <- colnames(d11)
hmp.ag.null <- rbind(d11,
                     d21)


d21 <- matrix(0, ncol = sum(ix.ag.other.tf), nrow = n.hmp)
d22 <- ag[, ix.ag.other]
colnames(d21) <- colnames(d22)
hmp.null.ag <- rbind(d21,
                     d22)

d.hmp.ag.full <- cbind(d.hmp.ag, hmp.ag.null, hmp.null.ag)
write.csv(d.hmp.ag.full, file = './_data/mm-HMP_data_082816-5-sites-level2-AG_L2_030317-with-bacteria-names-and-non-overlapping-cols.csv')





