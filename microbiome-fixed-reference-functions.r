#######################################################
# FUNCTIONS FOR THE HOTELLING'S T-TEST WITH BOOTSTRAP #
# #####################################################
# GHotelling
# in:
# data = matrix containg distances to one or both sites, so any of these rows (can be in any order):
# (x1, y1, NA, NA) = those subjects with samples from the first site only (group 1)
# (x1, y1, x2, y2) = those subjects with samples from both sites (group 3)
# (NA, NA, x2, y2) = those subjects with samples from the 2nd site only (group 2)
#
# nBoot - number of bootstrap resamples, sampling is within group
#
# print.details = do you want the intermediate calculations printed? (never printed in bootstrap loop)
#
# seed = random number generator seed, set before the bootstrap resampling
#
# out:
# returns the Hotelling t^2 statistic and bootstrap pvalue (if nBoot is not NULL)
#   accounting for correlation between samples from the same individual
GHotelling <- function(data, nBoot = NULL, print.details = F, seed = 1){

    colnames(data) <- c('x1', 'y1', 'x2', 'y2')
    n <- dim(data)[1]

    t2 <- GHotelling.helper(data, print.details = print.details)

    if(!is.null(nBoot)){

        mean.mx   <- matrix(rep(colMeans(data, na.rm = T), n), nrow = n, ncol = 4, byrow = T)
        data.null <- data - mean.mx

        boot.ix   <- get.boot.ix(data.null, nBoot, seed = seed)
        t2.boot   <- matrix(NA, nrow = nBoot, ncol = 1)

        for(i in 1:nBoot){
            t2.boot[i] <- GHotelling.helper(data = data.null[boot.ix[,i], ], print.details = F)
        }

        pval <- mean(t2.boot >= t2) # wrt original data
        if(print.details){
            print(paste('pval =', pval), q = F)
        }
    }else{
        pval <- NULL
    }
    return(list(t2 = t2, pval = pval))
}




# returns an list of 3 indices
# group1.ix = those subjects with samples from the first site only (x1, y1, NA, NA)
# group3.ix = those subjects with samples from both sites (x1, y1, x2, y2)
# group2.ix = those subjects with samples from the 2nd site only (NA, NA, x2, y2)
get.group.ix <- function(data){
    group.ix.a <- rowSums(!is.na(data[, 1:2]))
    group.ix.b <- rowSums(!is.na(data[, 3:4]))

    # the number of subjects in each group should be the same for group.ix and group.ix[boot.ix]
    group1.ix <- group.ix.a == 2 & group.ix.b == 0
    group3.ix <- group.ix.a == 2 & group.ix.b == 2
    group2.ix <- group.ix.a == 0 & group.ix.b == 2

    return(list(group1.ix = group1.ix,
                group3.ix = group3.ix,
                group2.ix = group2.ix))
}




# returns a matrix of indices for each bootstrap resample
# sampling is done within a group
get.boot.ix <- function(data.null, nBoot, seed = 1){
    set.seed(seed)

    rows <- 1:dim(data.null)[1]

    group.ix <- get.group.ix(data.null)

    # take out of list, since will be looping and accessing these
    group1.ix <- group.ix$group1.ix
    group3.ix <- group.ix$group3.ix
    group2.ix <- group.ix$group2.ix

    # the number of subjects in each group should be the same for group.ix and group.ix[boot.ix]
    n1 <- sum(group.ix$group1.ix)
    n3 <- sum(group.ix$group3.ix)
    n2 <- sum(group.ix$group2.ix)

    boot.ix <- matrix(NA, nrow = dim(data.null)[1], ncol = nBoot)
    for(i in 1:nBoot){
        boot.ix[,i] <- c(sample(rows[group1.ix], n1, replace = T),
                         sample(rows[group3.ix], n3, replace = T),
                         sample(rows[group2.ix], n2, replace = T))
    }
    return(boot.ix)
}




# returns the Hotelling t^2 statistic for a single dataset
GHotelling.helper <- function(data, print.details = F){ # nBoot not yet used

    group.ix <- get.group.ix(data)

    # the number of subjects in each group should be the same for group.ix and group.ix[boot.ix]
    n1 <- sum(group.ix$group1.ix)
    n3 <- sum(group.ix$group3.ix)
    n2 <- sum(group.ix$group2.ix)

    # covariance matrix
    a11 <- matrix(0, 2, 2)
    a12 <- matrix(0, 2, 2)
    a22 <- matrix(0, 2, 2)

    if(n1 > 0){
        a11[1,1] <- var(data[,'x1'], na.rm = T)
        a11[1,2] <- var(data[,'x1'], data[,'y1'], na.rm = T)
        a11[2,1] <- a11[1,2]
        a11[2,2] <- var(data[,'y1'], na.rm = T)
    }

    if(n2 > 0){
        a22[1,1] <- var(data[,'x2'], na.rm = T)
        a22[1,2] <- var(data[,'x2'], data[,'y2'], na.rm = T)
        a22[2,1] <- a22[1,2]
        a22[2,2] <- var(data[,'y2'], na.rm = T)
    }

    if(n3 > 0){
        a12[1,1] <- cov(x = data[,'x1'], y = data[,'x2'], use = 'pairwise.complete.obs')
        a12[1,2] <- cov(x = data[,'x1'], y = data[,'y2'], use = 'pairwise.complete.obs')
        a12[2,1] <- cov(x = data[,'y1'], y = data[,'x2'], use = 'pairwise.complete.obs')
        a12[2,2] <- cov(x = data[,'y1'], y = data[,'y2'], use = 'pairwise.complete.obs')
    }

    A <- matrix(0, nrow = 4, ncol = 4)
    A[1:2, 1:2] <- a11/(n1 + n3)
    A[3:4, 3:4] <- a22/(n3 + n2)
    A[1:2, 3:4] <- a12*(n3/((n1+n3)*(n3+n2)))
    A[3:4, 1:2] <- t(a12)*(n3/((n1+n3)*(n3+n2)))

    lambda <- matrix(c(1,0,-1,0,0,1,0,-1), byrow = F, ncol = 2)

    D <- matrix(c(mean(data[,'x1'], na.rm = T) - mean(data[,'x2'], na.rm = T),
                  mean(data[,'y1'], na.rm = T) - mean(data[,'y2'], na.rm = T)),
                ncol = 1)

    sigma <- t(lambda)%*%A%*%lambda

    t2 <- t(D)%*%solve(sigma)%*%D
    if(print.details){
        print('A, D, sigma, t^2:')
        print(A)
        print(D)
        print(sigma)
        print(t2)
    }

    return(as.numeric(t2))
}






##################################################
# FUNCTIONS FOR COMPUTING DISSIMILARITY MEASURES #
# ################################################

# new samples vs. (fixed reference)
# tax.level = taxonomy level, one of the following: "l2.phylum", "l3.class", "l4.order", "l5.family", "l6.genus", case does not matter
# d.new.filename = new data csv filename (together with path if needed) or rdata file. This is a matrix of relative abundances
            # in the new dataset (n samples x n leaves), with column names as in GreenGenes 13.8 (see GG_13_8_taxonomy.csv),
            # where n leaves = 2929 if leaf = genus, 1116 if leaf = family, 664 if leaf = order, 319 if leaf = class and 91 if
            # leaf = phyla. The matrix can contain other information in the columns, such as sample ID, etc.
            # The columns that are not relative abundance columns should be specified as an index passed to d.new.ix.col.not.rel.abu
# d.new.ix.col.not.rel.abu =  index of columns that are not relative abunance columns,
# measure = 'bc' or 'corr', default is bc
# automatically loads "taxonomy-GG13-8-and-HMP-references.rdata" from github, so needs a connection (taxonomy for GG13.8 file and HMP reference sets)
HMPdistance <- function(tax.level = NULL,                 # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                        d.new.filename = NULL,            # new data csv filename (together with path if needed to be read in by read.csv) or rdata file (with path)
                        d.new.ix.col.not.rel.abu = NULL,  # index of columns that are not relative abunance columns (ie. index of metadata)
                        measure = 'bc',                   # 'bc' for Bray-Curtis, 'corr' for 1-Pearson correlation
                        print.details = F){


    # these are the possible taxonomy levels, tax.level needs to be one of these
    tax.levels <- c('l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus')

    # load('taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata')
    source_data("https://github.com/NCI-biostats/microbiome-fixed-reference/blob/master/taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata?raw=True")
    # this should load gg.taxonomy.file, d.ref.stool, d.ref.nasal, d.ref.stool.info, d.ref.nasal.info
=
    # make sure things loaded correctly
    if(!length(c(grep('d.ref.stool', ls()), grep('d.ref.nasal', ls()), grep('gg.taxonomy.file', ls()))) >= 5){
        print('You need to be connected to the internet - file taxonomy-GG13-8-and-HMP-references.rdata needs to load from https://github.com/NCI-biostats/microbiome-fixed-reference. Exiting...', q = F)
        return(NULL)
    }
    d.gg.ra <- get.relative.abundance.matrix(tax.level = tax.level,                      # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                                             d.new.filename = d.new.filename,           # new data csv filename (together with path if needed) or rdata file
                                             d.new.ix.col.not.rel.abu = d.new.ix.col.not.rel.abu, # which columns in d.new are not relative abunance columns?
                                             gg.taxonomy.file = gg.taxonomy.file,
                                             d.ref.stool = d.ref.stool,
                                             d.ref.nasal = d.ref.nasal,
                                             tax.levels = tax.levels,
                                             print.details = print.details)


    dist.to.stool <- hmp.distance.helper(d.ref = d.ref.stool[[which(tax.levels == tax.level)]], # matrix of rel abu
                                         d.new = d.gg.ra,
                                         measure = measure)
    dist.to.nasal <- hmp.distance.helper(d.ref = d.ref.nasal[[which(tax.levels == tax.level)]], # matrix of rel abu
                                         d.new = d.gg.ra,
                                         measure = measure)

    hmp.dist <- cbind(dist.to.stool, dist.to.nasal)
    colnames(hmp.dist) <- c('dist.to.hmp.stool', 'dist.to.hmp.nasal')
    if(!is.null(rownames(d.gg.ra))){
        rownames(hmp.dist) <- rownames(d.gg.ra)
    }
    return(hmp.dist)
}






hmp.distance.helper <- function(d.ref, d.new, measure = 'bc'){
    n <- dim(d.new)[1]
    n.ref  <- dim(d.ref)[1]
    out <- matrix(NA, n, n.ref)

    if(!is.matrix(d.ref)){
        print('error in HMPdistance() - d.ref [matrix of relative abundances in the reference set] must be a matrix')
        return(NULL)
    }
    if(!is.matrix(d.new)){
        print('error in HMPdistance() - d.new [matrix of relative abundances for new samples] must be a matrix')
        return(NULL)
    }

    if(tolower(measure) == 'bc'){
        for(i in 1:n){
            for(j in 1:n.ref){
                out[i,j] <- get.bc(d.new[i,], d.ref[j,])
            }
        }
    }else if(tolower(measure) == 'corr' | tolower(measure) == 'cor'){
        for(i in 1:n){
            for(j in 1:n.ref){
                out[i,j] <- get.corr(d.new[i,], d.ref[j,], method = 'pearson')
            }
        }
    }else{
        print('error in HMPdistance() - option for measure not recognized, use either "bc" for Bray-Curtis or "corr" for Pearson correlation')
    }

    # return the dissimiliarity to the reference as a 1-column vector of length dim(d.new)[1]
    out.summary <- rowMeans(out)
    return(matrix(out.summary, ncol = 1))
}



# As in Bray-Curtis paper, does not assume that the vector x or y sums to 1
get.bc <- function(x, y){
    return(1 - sum(pmin.int(x, y)))
}



get.corr <- function(x, y, method = 'pearson'){
    return(1-cor(x, y, method = method))
}






###############################
# TAXONOMY MATCHING FUNCTIONS #
###############################

# modify the strings to make matching more universal/possible
# tax.org = vector of strings
cleanup.taxonomy.string <- function(tax.org){
    tax <- paste(tolower(tax.org), ';', sep = '') # append ; (needed for removing ;x__ at the lowest level)
    n <- length(tax)
    tax.clean <- vector('character', n)
    for(i in 1:n){
        curr1 <- gsub('other;{1}', ';' , tax[i], ignore.case = T)        # get rid of Other or other
        curr2 <- gsub('[[:alpha:]]__;{1}', ';', curr1, ignore.case = T) # get rid of empty fields such as p__ or o__ etc., leave the spacing in, so ;p__; becomes ;;
        curr3 <- gsub('i{2}', 'i', curr2, ignore.case = T) # ii -> i
        curr4 <- gsub('\\[', '', curr3)
        curr5 <- gsub('\\]', '', curr4)
        tax.clean[i] <- gsub('Root;', '', curr5, ignore.case = T)       # get rid of Root;
    }
    if(substr(tax.clean[1], 1, 3) != 'k__'){ # check the first, all the rest will be the same
        ix <- grepl('archaeota', tax.clean) > 0
        tax.clean[ix]  <- paste('k__archaea;', tax.clean[ix], sep = '') # if no Kingdom specified, assume it is 'k__Bacteria' unless p__Euryarchaeota
        tax.clean[!ix] <- paste('k__bacteria;', tax.clean[!ix], sep = '')
    }
    return(tax.clean)
}




classify.not.found.as.missing <- function(tax.level, d.tax, gg.tax,
                                          tax.levels = c('l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'),
                                          tax.prefixes = c('p__', 'c__', 'o__', 'f__', 'g__'), print.details = F){
    curr.level <- which(tax.level == tax.levels)
    ix.missing <- is.na(match(x = d.tax, table = gg.tax))

    d.tax.out <- d.tax
    while(sum(ix.missing) > 0 & curr.level > 0){
        d.tax.missing <- d.tax.out[ix.missing]
        d.tax.updated <- gsub(paste(tax.prefixes[curr.level], '[[:alnum:]]*-*[.]*[[:alnum:]]*-*[.]*[[:alnum:]]*;', sep = ''), ';', d.tax.missing, ignore.case = T)
        if(print.details){
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            print(paste('in classify.not.found.as.missing(), starting with ', tax.level, ', currently at level ', tax.levels[curr.level], sep = ''), q = F)
            print('starting taxonomy for non-matched strings                                  changed to (by removing the last level', q = F)
            print(cbind(d.tax.out[ix.missing], d.tax.updated))
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        }
        d.tax.out[ix.missing] <- d.tax.updated
        ix.missing <- is.na(match(x = d.tax.out, table = gg.tax))
        curr.level <- curr.level - 1
    }
    return(d.tax.out)
}





get.gg.matched.ra.matrix <- function(gg.tax, d.tax.fixed, d.new, d.new.ix.col.not.rel.abu, gg.tax.org)
{
    ix.match <- match(x = d.tax.fixed, table = gg.tax)

    d.new.ra <- d.new[, -d.new.ix.col.not.rel.abu] # matrix of rel. abu from the user

    d.out <- matrix(0, nrow = dim(d.new)[1], ncol = length(gg.tax))
    colnames(d.out) <- gg.tax.org
    if(length(ix.match) != dim(d.new.ra)[2]){print('in get.gg.matched.ra.matrix() - there is something wrong')}
    for(i in 1:length(ix.match)){
        d.out[,ix.match[i]] <- d.out[,ix.match[i]] + d.new.ra[, i]
    }
    return(d.out)
}





get.relative.abundance.matrix <- function(tax.level,                # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                                          d.new.filename,           # new data csv filename (together with path if needed) or rdata file
                                          d.new.ix.col.not.rel.abu, # which columns in d.new are not relative abunance columns?
                                          gg.taxonomy.file = gg.taxonomy.file,
                                          d.ref.stool = d.ref.stool,
                                          d.ref.nasal = d.ref.nasal,
                                          tax.levels,               # = c('l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'),
                                          print.details)            # = F) # order of columns corresponds to columns in gg.taxonomy.file
{
    # GreanGenes taxonomy
    if(sum(tolower(tax.level) %in% tax.levels) != 1){
        print("error in match.columns() - tax.level should be one of: 'L2.phylum', 'L3.class', 'L4.order', 'L5.family', 'L6.genus. Case does not matter.'", q = F)
    }
    gg.tax.org <- levels(gg.taxonomy.file[[which(tax.levels == tax.level)]]) # extract the list of taxonomy names at the correct level
    if(gg.tax.org[1] == ""){ # remove the empty level
        gg.tax.org <- gg.tax.org[-1]
    }
    gg.tax <- cleanup.taxonomy.string(gg.tax.org)         # returns a cleaned up taxonomy


    # taxonomy used in the new data
    # load or read in the data
    if(substr(d.new.filename, nchar(d.new.filename)-2, nchar(d.new.filename)) == 'csv'){
        d.new     <- read.csv(d.new.filename, check.names = F)
    }else{
        load(d.new.filename)
    }

    # check that the files read in correctly, or were provided in a correct format.
    if(sum(grepl(';', names(d.new))) < 1){
        print(paste('Make sure the names are formatted according to GreenGenes 13.8. See GG_13_8_taxonomy.csv file in Github under ',
                    github.path, '. Specifically with ; as a delimiter. If using read.csv, use option check.names = F.', sep = ''))
    }

    if(is.null(d.new.ix.col.not.rel.abu)){
        d.tax.org <- colnames(d.new)
    }else{
        d.tax.org <- colnames(d.new)[-d.new.ix.col.not.rel.abu] # remove the columns that do not correspond to the relative abundance matrix
    }
    d.tax     <- cleanup.taxonomy.string(d.tax.org)

    # many of the HMP taxonomic names are not found in GG 13.8 taxonomy
    # for example, for phyla: "Root.p__CCM11b"    "Root.p__Thermi"        "Root.p__WPS.2"         "Root.p__ZB2"
    ix.missing <- is.na(match(x = d.tax, table = gg.tax))

    if(print.details){
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print(paste(sum(ix.missing), ' taxonomies in the dataset not found in the GG 13.8 taxonomy at level ', tax.level,
                    '. They will be classified as missing at the lowest possible level. They are:', sep = ''), q = F)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print(d.tax[which(ix.missing)])
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    }

    d.tax.fixed <- classify.not.found.as.missing(tax.level, d.tax, gg.tax, print.details = print.details)

    d.new.gg <- get.gg.matched.ra.matrix(gg.tax, d.tax.fixed, d.new, d.new.ix.col.not.rel.abu, gg.tax.org)
    if(!is.null(rownames(d.new))){
        rownames(d.new.gg) <- rownames(d.new)
    }
    return(d.new.gg)
}




