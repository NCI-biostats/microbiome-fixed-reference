##################################################
# FUNCTIONS FOR COMPUTING DISSIMILARITY MEASURES #
# ################################################

# new samples vs. (fixed reference)
# tax.level = taxonomy level, one of the following: "l2.phylum", "l3.class", "l4.order", "l5.family", "l6.genus", case does not matter
# d.new.filename = new data csv filename (together with path if needed) or rdata file.
#                  This is either a matrix of relative abundances in the new dataset (n samples x n phyla), with columnames as in GreenGenes 13.8 (see GG_13_8_taxonomy.csv)
#                  or it can contain other information in the first k columns (d.new.first.k.col.not.rel.abu), those need to be specified and will be preserved.
# measure = 'bc' or 'corr', default is bc

HMPdistance <- function(tax.level = NULL,                # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                        d.new.filename = NULL,            # new data csv filename (together with path if needed) or rdata file
                        d.new.ix.col.not.rel.abu = NULL, # index of columns that are not relative abunance columns?
                        measure = 'bc',
                        tax.levels = c('l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'),

                        github.path = '',                # with / at the end
                        gg.taxonomy.file = 'GG_13_8_taxonomy.csv',
                        reference.sets.filename = 'HMP_L2_thru_L6_030317_GG13_8_ref_sets_stool_nasal.rdata',   # d.ref.l2.phylum.nasal, d.ref.l2.phylum.stool, d.ref.l3.class.nasal, d.ref.l3.class.stool,
                                                                                                        # d.ref.l4.order.nasal, d.ref.l4.order.stool, d.ref.l5.family.nasal, d.ref.l5.family.stool,
                                                                                                        # d.ref.l6.genus.nasal, d.ref.l6.genus.stool, d.ref.nasal.info, d.ref.stool.info
                        print.details = F){




    d.gg.ra <- get.relative.abundance.matrix(tax.level = tax.level,                      # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                                              d.new.filename = d.new.filename,           # new data csv filename (together with path if needed) or rdata file
                                              d.new.ix.col.not.rel.abu = d.new.ix.col.not.rel.abu, # which columns in d.new are not relative abunance columns?
                                              github.path = github.path,                 # with / at the end
                                              gg.taxonomy.file = gg.taxonomy.file,
                                              tax.levels = tax.levels,
                                              print.details = print.details)


    load(paste(github.path, reference.sets.filename, sep = ''))  # "d.ref.nasal"      "d.ref.nasal.info" "d.ref.stool"      "d.ref.stool.info"


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



get.relative.abundance.matrix <- function(tax.level,# = NULL,                # one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
                                          d.new.filename,# = NULL,           # new data csv filename (together with path if needed) or rdata file
                                          d.new.ix.col.not.rel.abu,# = NULL, # which columns in d.new are not relative abunance columns?
                                          github.path,#  = '',                # with / at the end
                                          gg.taxonomy.file,#  = 'GG_13_8_taxonomy.csv',
                                          tax.levels, # = c('l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'),
                                          print.details)# = F) # order of columns corresponds to columns in gg.taxonomy.file
{
    tax.level <- tolower(tax.level)

    # GreanGenes taxonomy
    gg.tax.all <- read.csv(paste(github.path, gg.taxonomy.file, sep = '')) # read in as factors
    if(sum(tolower(tax.level) %in% tax.levels) != 1){
        print("error in match.columns() - tax.level should be one of: 'L2.phylum', 'L3.class', 'L4.order', 'L5.family', 'L6.genus. Case does not matter.'", q = F)
    }
    gg.tax.org <- levels(gg.tax.all[[which(tax.levels == tax.level)]]) # extract the list of taxonomy names at the correct level
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





















#######################################################
# FUNCTIONS FOR THE HOTELLING'S T-TEST WITH BOOTSTRAP #
# #####################################################

#
# hotelling.test(x = d.a, y = d.b, shrinkage = shrinkage, perm = perm, B = B)
#
# HotellingsT2(X = d.a, Y = NULL, mu = mu, test = test)
#
#
# hotelling.test.fixed.ref <- function(x, y, )



hotelling.test.fixed.reference <- function()


get.row.hotelling.t2.and.boot <- function(d, a.site, a.ref = T, b.site, b.ref = F, measure, B, main.paper = F){
    if(measure == 'bc'){
        x <- with(d, cbind(bc1[ref == a.ref & site == a.site], bc2[ref == a.ref & site == a.site]))
        y <- with(d, cbind(bc1[ref == b.ref & site == b.site], bc2[ref == b.ref & site == b.site]))
    }else if(measure == 'corr.p'){
        x <- with(d, cbind(corr.p1[ref == a.ref & site == a.site], corr.p2[ref == a.ref & site == a.site]))
        y <- with(d, cbind(corr.p1[ref == b.ref & site == b.site], corr.p2[ref == b.ref & site == b.site]))
    }else if(measure == 'corr.s'){
        x <- with(d, cbind(corr.s1[ref == a.ref & site == a.site], corr.s2[ref == a.ref & site == a.site]))
        y <- with(d, cbind(corr.s1[ref == b.ref & site == b.site], corr.s2[ref == b.ref & site == b.site]))
    }

    out <- hotelling.test(x, y)

    data.out <- get.boot.dataset(d, site1 = a.site, site2 = b.site, site.1.ref = a.ref, site.2.ref = b.ref, measure)
    boot.out <- get.boot.t2.p.val(data.org = data.out$data.org, data.out$group.ix, B = B)
    if(main.paper){
        row.out <- matrix(c(sum(data.out$group.ix == 1), sum(data.out$group.ix == 3), sum(data.out$group.ix == 2), out$stats$statistic, out$pval, boot.out$statistic, boot.out$pval), nrow = 1)
    }else{
        row.out <- matrix(c(out$stats$statistic, out$pval, boot.out$statistic, boot.out$pval, sum(data.out$group.ix == 1), sum(data.out$group.ix == 3), sum(data.out$group.ix == 2)), nrow = 1)
    }
    rownames(row.out) <- paste(a.site, ' vs ', b.site, sep = '')
    return(row.out)


}





# site1.set = 'ref' or 'test'
get.boot.dataset <- function(d, site1 = 'skin', site2 = 'vaginal', site.1.ref = T, site.2.ref = T, measure){

    attach(d)
    site1.ix <- site == site1 & ref == site.1.ref
    site2.ix <- site == site2 & ref == site.2.ref

    site1.id <- id[site1.ix]
    site2.id <- id[site2.ix]

    u.ids <- matrix(unique(c(site1.id, site2.id)), ncol = 1)

    u.ids.group <- matrix(NA, ncol = 1, nrow = length(u.ids))
    for(i in 1:length(u.ids)){
        ix1 <- u.ids[i] %in% site1.id
        ix2 <- u.ids[i] %in% site2.id
        if(sum(ix1) == 1 & sum(ix2) == 1){
            u.ids.group[i] <- 3
        }else if (sum(ix1) == 1 & sum(ix2) == 0){
            u.ids.group[i] <- 1
        }else if (sum(ix1) == 0 & sum(ix2) == 1){
            u.ids.group[i] <- 2
        }
    }

    group1.ix  <- (id %in% u.ids[u.ids.group == 1]) & site1.ix
    group2.ix  <- (id %in% u.ids[u.ids.group == 2]) & site2.ix
    group31.ix <- (id %in% u.ids[u.ids.group == 3]) & site1.ix
    group32.ix <- (id %in% u.ids[u.ids.group == 3]) & site2.ix

    detach(d)

    n1 <- sum(group1.ix)
    n2 <- sum(group2.ix)
    n3 <- sum(group31.ix)

    group.ix <- matrix(c(rep(1, n1),
                         rep(3, n3),
                         rep(2, n2)), ncol = 1)

    measure.col.ix <- grep(measure, names(d))

    data.group1  <- as.matrix(d[group1.ix, measure.col.ix])
    data.group2  <- as.matrix(d[group2.ix, measure.col.ix])
    data.group31 <- as.matrix(d[group31.ix, measure.col.ix])
    data.group32 <- as.matrix(d[group32.ix, measure.col.ix])

    if(length(measure.col.ix) == 1){ # one references
        data.org <- matrix(NA, ncol = 2, nrow = length(u.ids.group))
        if(n1 > 0){
            data.org[1:n1, 1] <- data.group1
        }
        if(n3 > 0){
            data.org[(n1+1):(n1+n3), 1] <- data.group31
            data.org[(n1+1):(n1+n3), 2] <- data.group32
        }
        if(n2 > 0){
            data.org[(n1+n3+1):(n1+n3+n2), 2] <- data.group2
        }
        colnames(data.org) <- c('x1', 'x2')
    }else if(length(measure.col.ix) == 2){ # two references
        data.org <- matrix(NA, ncol = 4, nrow = length(u.ids.group))
        if(n1 > 0){
            data.org[1:n1, 1:2] <- data.group1
        }
        if(n3 > 0){
            data.org[(n1+1):(n1+n3), 1:2] <- data.group31
            data.org[(n1+1):(n1+n3), 3:4] <- data.group32
        }
        if(n2 > 0){
            data.org[(n1+n3+1):(n1+n3+n2), 3:4] <- data.group2
        }
        colnames(data.org) <- c('x1', 'y1', 'x2', 'y2')
    }
    return(list(data.org = data.org, group.ix = group.ix))
}




# data = matrix x1, y1, x2, y2: x = ref 1, y = ref 2, 1 = site 1, 2 site 2
# group.ix = 1,2,3: 1 = samples from site 1 only, 2 = site 2 only, 3 = site 1 and 2
get.boot.t2.p.val.one.ref <- function(data.org, group.ix, B){
    x1.mn <- mean(data.org[,1], na.rm = T) # mean(x1, na.rm = T)
    x2.mn <- mean(data.org[,2], na.rm = T) # mean(x2, na.rm = T)

    x.mn <- mean(c(data.org[,1], data.org[,2]), na.rm = T) # mean(c(x1, x2), na.rm = T)

    mean.mx <- matrix(rep(c(x1.mn, x2.mn), dim(data.org)[1]), nrow = dim(data.org)[1], ncol = 2, byrow = T)

    data.null.mx <- data.org - mean.mx

    colnames(data.null.mx) <- c('x1', 'x2')

    t2.obs.org    <- get.boot.t2(data.org, group.ix = group.ix, boot.ix = NULL, one.ref = T)
    t2.obs.org.mx <- matrix(t2.obs.org, nrow = B, ncol = 1)

    t2.boot <- matrix(NA, nrow = B, ncol = 1)

    for(i in 1:B){
        boot.ix    <- get.boot.ix(group.ix)
        t2.boot[i] <- get.boot.t2(data = data.null.mx, group.ix = group.ix, boot.ix = boot.ix, one.ref = T)
    }

    boot.p.val <- mean(t2.boot >= t2.obs.org.mx) # wrt original data
    return(list(statistic = t2.obs.org, pval = boot.p.val))
}



# data = matrix x1, y1, x2, y2: x = ref 1, y = ref 2, 1 = site 1, 2 site 2
# group.ix = 1,2,3: 1 = samples from site 1 only, 2 = site 2 only, 3 = site 1 and 2
get.boot.t2.p.val.two.ref <- function(data.org, group.ix, B){

    x1.mn <- mean(data.org[,1], na.rm = T) # mean(x1, na.rm = T)
    y1.mn <- mean(data.org[,2], na.rm = T) # mean(y1, na.rm = T)
    x2.mn <- mean(data.org[,3], na.rm = T) # mean(x2, na.rm = T)
    y2.mn <- mean(data.org[,4], na.rm = T) # mean(y2, na.rm = T)

    x.mn <- mean(c(data.org[,1], data.org[,3]), na.rm = T) # mean(c(x1, x2), na.rm = T)
    y.mn <- mean(c(data.org[,2], data.org[,4]), na.rm = T) # mean(c(y1, y2), na.rm = T)


    mean.mx <- matrix(rep(c(x1.mn, y1.mn, x2.mn, y2.mn), dim(data.org)[1]), nrow =  dim(data.org)[1], ncol = 4, byrow = T)

    data.null.mx <- data.org - mean.mx

    colnames(data.null.mx) <- c('x1', 'y1', 'x2', 'y2')

    t2.obs.org    <- get.boot.t2(data.org, group.ix = group.ix, boot.ix = NULL)
    t2.obs.org.mx <- matrix(t2.obs.org, nrow = B, ncol = 1)

    t2.boot <- matrix(NA, nrow = B, ncol = 1)

    for(i in 1:B){
        boot.ix    <- get.boot.ix(group.ix)
        t2.boot[i] <- get.boot.t2(data = data.null.mx, group.ix = group.ix, boot.ix = boot.ix)
    }

    boot.p.val <- mean(t2.boot >= t2.obs.org.mx) # wrt original data
    return(list(statistic = t2.obs.org, pval = boot.p.val))
}
