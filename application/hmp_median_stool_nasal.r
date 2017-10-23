rm(list = ls())
#setwd('change this if needed')
source('code-microbiome-helper-functions.r')

d.full <- read.csv('data/yw-HMP_data_082816-5-sites-level2-mm.csv')

# list of unique id's
u.id <- unique(d.full$id)
n.id <- length(u.id)

set.seed(1)

#split samples in two groups (ref set and test set)
s1    <- sample(1:n.id, size = n.id/2, replace = F)
s1.tf <- rep(F, n.id) # T/F index
s1.tf[s1] <- T

# list of unique train and test id's
id.ref  <- u.id[s1.tf]   # ref set, 97 individuals
id.test <- u.id[!s1.tf]  # test set, 98 individuals

attach(d.full)
d <- data.frame(id = id,
                sid = sample.id,
                ref = id %in% id.ref,
                female = sex == 'female',
                site = site,
                stool = site == 'stool',
                saliva = site == 'saliva',
                skin = site == 'skin',
                nasal = site == 'nasal',
                vaginal = site == 'vaginal')


d.ref.ix1 = d$ref & d$stool # this is the index to the rows of the reference stool samples
d.ref.ix2 = d$ref & d$nasal # this is the index to the rows of the reference nasal samples 

dp <- data.matrix(d.full[, -(1:4)]) # remove sample.id, id, sex, site, leaving only the relative abundances.
detach(d.full)

out <- table(!d$ref, d$site)
rownames(out) <- c('Set (A) Ref', 'Set (B)')

out

# nasal saliva skin stool vaginal
# Set (A) Ref    74     75   56    92      33
# Set (B)        72     83   62    87      47

#using median distance. 
do <- d
d$bc     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1)
d$corr.s <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1) # spearman
d$corr.p <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1) # pearson

#save(d, file = 'mm-HMP_data_082816-5-sites-level2-stool-ref-half-sample-median.rdata')
rm(d)

d <- do
d$bc     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1)
d$corr.s <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1) # spearman
d$corr.p <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1) # pearson

#save(d, file = 'mm-HMP_data_082816-5-sites-level2-nasal-ref-half-sample-median.rdata')
rm(d)

d <- do
d$bc1     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1)
d$corr.s1 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1) # spearman
d$corr.p1 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix1, stat = 'median'), ncol = 1) # pearson
d$corr.s2 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1) # spearman
d$bc2     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1)
d$corr.p2 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix2, stat = 'median'), ncol = 1) # pearson

#save(d, file = 'mm-HMP_data_082816-5-sites-level2-stool-nasal-refs-half-sample-median.rdata')
rm(d)