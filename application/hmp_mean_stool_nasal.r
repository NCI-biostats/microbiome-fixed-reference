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


do <- d

#distance with stool references.  
#d$bc     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix1), ncol = 1)
#d$corr.s <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix1), ncol = 1) # spearman
#d$corr.p <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix1), ncol = 1) # pearson

#mn.stool.ref.bc <- get.ref.means(d, dp, measure = 'bc', ref.ix = d.ref.ix1)
#mn.stool.ref.sc <- get.ref.means(d, dp, measure = 'corr.s', ref.ix = d.ref.ix1)
#mn.stool.ref.pc <- get.ref.means(d, dp, measure = 'corr.p', ref.ix = d.ref.ix1)

#save(d, mn.stool.ref.bc, mn.stool.ref.sc, mn.stool.ref.pc, file = 'mm-HMP_data_082816-5-sites-level2-stool-ref-half-sample-mean.rdata')
#rm(d)

#d <- do

#distances with nasal references
#d$bc  <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix2), ncol = 1)
#d$corr.s <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix2), ncol = 1) # spearman
#d$corr.p <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix2), ncol = 1) # pearson

#mn.nasal.ref.bc <- get.ref.means(d, dp, measure = 'bc', ref.ix = d.ref.ix2)
#mn.nasal.ref.sc <- get.ref.means(d, dp, measure = 'corr.s', ref.ix = d.ref.ix2)
#mn.nasal.ref.pc <- get.ref.means(d, dp, measure = 'corr.p', ref.ix = d.ref.ix2)

#save(d, mn.nasal.ref.bc, mn.nasal.ref.sc, mn.nasal.ref.pc, file = 'mm-HMP_data_082816-5-sites-level2-nasal-ref-half-sample-mean.rdata')
#rm(d)

#(YW note: why do it again? May no need for code before)
d <- do
d$bc1     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix1), ncol = 1)
d$corr.s1 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix1), ncol = 1) # spearman
d$corr.p1 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix1), ncol = 1) # pearson
d$bc2     <- as.matrix(get.distances(d, dp, measure = 'bc.mg', ref.ix = d.ref.ix2), ncol = 1)
d$corr.s2 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'spearman', ref.ix = d.ref.ix2), ncol = 1) # spearman
d$corr.p2 <- as.matrix(get.distances(d, dp, measure = 'corr', method = 'pearson', ref.ix = d.ref.ix2), ncol = 1) # pearson

# means calculation (BC and Pearson) for the refence set for each site wrt stool and nasal references
# (diag removed from mean calculation for stool and nasal)
sites <- c('stool', 'nasal', 'skin', 'saliva', 'vaginal')

# ~~~~~~~~~~~~~~~
# get the means for each site with respect to stool and nasal for bray curtis distance
mu.bc <- matrix(NA, ncol = 2, nrow = length(sites))  # columns: wrt stool, wrt nasal
rownames(mu.bc) <- sites
colnames(mu.bc) <- c('ref.stool', 'ref.nasal')
for(i in 1:5){
    ix <- d$ref & d$site == sites[i]
    x1.bc <- d$bc1[ix]      # wrt stool
    mu.bc[i,1] <- mean(x1.bc)
    x2.bc <- d$bc2[ix]      # wrt nasal
    mu.bc[i,2] <- mean(x2.bc)
}
# fix up the means of stool wrt stool and nasal wrt nasal by removing the diagonal zeros.
mu.bc[1,1] <- get.ref.stats(d, dp, measure = 'bc', ref.ix = d.ref.ix1, stat = 'mean') # stool to stool
mu.bc[2,2] <- get.ref.stats(d, dp, measure = 'bc', ref.ix = d.ref.ix2, stat = 'mean') # nasal to nasal

mu.bc
# ~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~
# get the means for each site with respect to stool and nasal for pearson distance
mu.pc <- matrix(NA, ncol = 2, nrow = length(sites)) # columns: wrt stool, wrt nasal
rownames(mu.pc) <- sites
colnames(mu.pc) <- c('ref.stool', 'ref.nasal')
for(i in 1:5){
    ix <- d$ref & d$site == sites[i]
    x1.pc <- d$corr.p1[ix]      # wrt stool
    mu.pc[i,1] <- mean(x1.pc)
    x2.pc <- d$corr.p2[ix]      # wrt nasal
    mu.pc[i,2] <- mean(x2.pc)
}
# fix up the means of stool wrt stool and nasal wrt nasal by removing the diagonal zeros.
mu.pc[1,1] <- get.ref.stats(d, dp, measure = 'corr.p', ref.ix = d.ref.ix1, stat = 'mean') # stool to stool
mu.pc[2,2] <- get.ref.stats(d, dp, measure = 'corr.p', ref.ix = d.ref.ix2, stat = 'mean') # nasal to nasal

mu.pc
# ~~~~~~~~~~~~~~~

#save(d, mu.pc, mu.bc, file = 'mm-HMP_data_082816-5-sites-level2-stool-nasal-refs-half-sample-mean.rdata')
#save(d, dp, mu.pc, mu.bc, file = 'mm-HMP_data_082816-5-sites-level2-stool-nasal-refs-half-sample-mean-raw-rel-abu.rdata')

rm(d)

