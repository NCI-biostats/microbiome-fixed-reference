
# As in Bray-Curtis paper
get.bc <- function(x, y){
  mins <- apply(cbind(x, y), 1, min)
  return(1 - 2*sum(mins)/(sum(x) + sum(y)))
}

get.bc.qiime <- function(x, y){
  num <- sum(abs(x-y))
  den <- sum(x+y)
  return(num/den)
}


get.bc.mg <- function(x, y){
  mins <- apply(cbind(x, y), 1, min)
  return(1 - sum(mins))
}


get.corr <- function(x, y, method = 'spearman'){
  return(1-cor(x, y, method = method))
}


get.covariance <- function(x, y){
  return(1-cov(x, y))
}



# all samples vs. (fixed reference)
# d = dataset    (n samples x (id, sid, ref, female, site, stool, saliva, skin, nasal, vaginal, bc, corr.s, corr.p) - see below)
# dp = matrix of relative abundances (n samples x n phyla)
# measure = 'bc.mg', 'corr.p', 'corr.s'
# method = 'spearman' or 'pearson'
# ref.ix = reference set index, in our analysis this is either: set A (ref) & stool, or set A and nasal
# stat = 'mean' or 'median', others can be added
# return.full.mx = return the full matrix of distances (n samples x n samples in ref.ix),
# if return.full.mx = F, returns the stat for each row (n samples x 1)
get.distances <- function(d, dp, measure = 'bc.mg', method = 'spearman', ref.ix, stat = 'mean', return.full.mx = F){
  n <- dim(d)[1]

  dp.ref <- dp[ref.ix, ]
  n.ref  <- dim(dp.ref)[1]
  out <- matrix(NA, n, n.ref)

  if(measure == 'bc.mg'){
    for(i in 1:n){
      for(j in 1:n.ref){
        out[i,j] <- get.bc.mg(dp[i,], dp.ref[j,])
      }
    }
  }else if(measure == 'corr'){
    for(i in 1:n){
      for(j in 1:n.ref){
        out[i,j] <- get.corr(dp[i,], dp.ref[j,], method = method)
      }
    }
  }
  if(return.full.mx){
      return(out)
  }

  if(stat == 'mean'){
      return(rowMeans(out))
  }else if(stat == 'median'){
      return(apply(out, 1, median))
  }
}




get.distances.wrt.ref.mean.vec <- function(d, dp, measure = 'bc.mg', method = 'spearman', ref.ix, stat = NULL, return.full.mx = F){
    n <- dim(d)[1]

    dp.ref.mx     <- dp[ref.ix, ]
    dp.ref.mn.vec <- apply(dp.ref.mx, 2, mean)
    out           <- rep(NA, n)
    if(measure == 'bc.mg'){
        for(i in 1:n){
            out[i] <- get.bc.mg(dp[i,], dp.ref.mn.vec)
        }
    }else if(measure == 'corr'){
        for(i in 1:n){
            out[i] <- get.corr(dp[i,], dp.ref.mn.vec, method = method)
        }
    }
    return(out)

}





t.test.wrap <- function(d, a.site, a.ref, b.site, b.ref, measure){
  d.a <- d[d$ref == a.ref & d$site == a.site, measure]
  d.b <- d[d$ref == b.ref & d$site == b.site, measure]
  t.out <- t.test(d.a, d.b)
  return(c(t.out$statistic, t.out$p.val))
}





get.row.hotelling <- function(d, a.site, a.ref = F, b.site, b.ref = T, measure, B = 1000, perm = F, shrinkage = F){
    d.a <- d[d$ref == a.ref & d$site == a.site, paste(measure, 1:2, sep = '')]
    d.b <- d[d$ref == b.ref & d$site == b.site, paste(measure, 1:2, sep = '')]

    out <- hotelling.test(x = d.a, y = d.b, shrinkage = shrinkage, perm = perm, B = B)
    return(c(out$stats$statistic, out$pval))
}




get.row.hotelling.one.sample <- function(d, a.site, a.ref = F, mu, measure, test = 'chi'){
    d.a <- d[d$ref == a.ref & d$site == a.site, paste(measure, 1:2, sep = '')]

    out <- HotellingsT2(X = d.a, Y = NULL, mu = mu, test = test)
    return(c(out$statistic, out$p.value))
}



# d = full dataset, including distances
# a.site, b.site = text, sites for a and b
# a.ref, b.ref = T/F
# measure = text, either 'bc' or 'corr'
t.test.jack <- function(d, a.site, a.ref, b.site, b.ref, measure){
  # sort d by id
  ix.sort <- sort(d$id, index.return = T)$ix
  d.s <- d[ix.sort, ]

  d.curr <- data.frame(id = unique(d.s$id), a = NA, b = NA)
  d.a <- d.s[d.s$ref == a.ref & d.s$site == a.site, c('id', measure)]
  d.b <- d.s[d.s$ref == b.ref & d.s$site == b.site, c('id', measure)]

  ix.a <- d.curr$id %in% d.a$id # this should be correct since the id's are sorted
  ix.b <- d.curr$id %in% d.b$id

  d.curr$a[ix.a] <- d.a[, measure]
  d.curr$b[ix.b] <- d.b[, measure]

  # now remove the subjects with no data in either sample
  ix.na <- is.na(d.curr$a) & is.na(d.curr$b)
  d.curr <- d.curr[!ix.na, ]
  n <- dim(d.curr)[1]

  # ready to do the jackknife
  diffs <- matrix(NA, n)
  for(i in 1:n){
    d.i <- d.curr[-i,]
    diffs[i] <- mean(d.i$a, na.rm = T) - mean(d.i$b, na.rm = T)
  }

  var.jack <- (n-1)/n*sum((diffs - mean(diffs))^2)
  if(var.jack == 0){
      t.stat <- 0
  }else{
      t.stat <- mean(diffs)/sqrt(var.jack)
  }
  if(t.stat < 0){
    p.val <- 2*pt(t.stat, df = n - 2, lower.tail = T)
  }else{
    p.val <- 2*pt(t.stat, df = n - 2, lower.tail = F)
  }
  return(c(t.stat, p.val))
}


t.test.wrap.one.sample <- function(d, a.site, a.ref, mu, measure, square.stat = F){
  d.a <- d[d$ref == a.ref & d$site == a.site, measure]
  t.out <- t.test(x = d.a, y = NULL, mu = mu)
  if(square.stat){
      return(c(t.out$statistic^2, t.out$p.val))
  }else{
      return(c(t.out$statistic, t.out$p.val))
  }
}


get.row.one.sample <- function(d, a.site, a.ref, mu, measure){
  d.a <- d[d$ref == a.ref & d$site == a.site, measure]
  out <- matrix(c(get.mean.median.iqr(d.a), # mean, median, IQR
                  mu,
                  t.test.wrap.one.sample(d, a.site, a.ref, mu, measure)), nrow = 1)
  rownames(out) <- a.site
  return(out)
}



get.row.no.boot <- function(d, a.site, a.ref, b.site, b.ref, measure){
    d.a <- d[d$ref == a.ref & d$site == a.site, measure]
    d.b <- d[d$ref == b.ref & d$site == b.site, measure]
    out <- matrix(c(summary(d.a)[c(4,3,2,5)], # mean, median, IQR
                    summary(d.b)[c(4,3,2,5)],
                    t.test.wrap(d, a.site, a.ref, b.site, b.ref, measure),
                    t.test.jack(d, a.site, a.ref, b.site, b.ref, measure)), nrow = 1)
    rownames(out) <- paste(a.site, ' vs ', b.site, sep = '')
    return(out)
}


get.mean.median.iqr <- function(x){
    return(c(mean(x),
             quantile(x, c(.5, .25, .75))))
}


get.mean.median <- function(x){
    return(c(mean(x), median(x)))
}


get.mean.sd <- function(x){
    return(c(mean(x, na.rm = T), sd(x, na.rm = T)))
}






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








get.row.table2.one.sample.bootstrap <- function(d, a.site, a.ref = F, ref.mu, measure, B = 10000){

    x <- d[d$ref == a.ref & d$site == a.site, grep(measure, names(d))] # 2 columns

    n <- dim(x)[1] # n subjects in test set for a.site

    sigma.hat <- var(x) # 2x2

    x.org.mn <- matrix(apply(x, 2, mean), nrow = 1) # means of the original, 1x2

    d.bar <- matrix(x.org.mn - ref.mu, nrow = 1) # 1x2 matrix

    t2.11 <- n * d.bar[1]^2/sigma.hat[1,1]
    t2.22 <- n * d.bar[2]^2/sigma.hat[2,2]
    t2    <- n * d.bar%*%solve(sigma.hat)%*%t(d.bar)


    # bootstrap
    t2.boot.11 <- rep(NA, B)
    t2.boot.22 <- rep(NA, B)
    t2.boot    <- rep(NA, B)

    # x.org.mn made into a matrix so they are correctly subtracted from x.boot
    x.org.mn.mx <- matrix(rep(x.org.mn, n), nrow = n, byrow = T)
    for(i in 1:B){
        x.boot <- as.matrix(x[sample(1:n, size = n, replace = T), ]) #  n x 2

        d.star <- x.boot - x.org.mn.mx  # n x 2

        sigma.hat.boot <- var(d.star)   # 2 x 2

        d.bar.boot <- matrix(apply(d.star, 2, mean), nrow = 1) # 1 x 2 matrix

        t2.boot.11[i]   <- n * d.bar.boot[1]^2/sigma.hat.boot[1,1]
        t2.boot.22[i]   <- n * d.bar.boot[2]^2/sigma.hat.boot[2,2]
        t2.boot[i]      <- n * d.bar.boot%*%solve(sigma.hat.boot)%*%t(d.bar.boot)
    }

    p.val.11 <- mean(t2.boot.11 > t2.11)
    p.val.22 <- mean(t2.boot.22 > t2.22)
    p.val    <- mean(t2.boot    > as.numeric(t2))

    return(c(n, t2.11, p.val.11, t2.22, p.val.22, t2, p.val))
}




get.row.table3 <- function(d, a.site, a.ref = F, b.site, b.ref = T, measure, B = 1000, one.ref = F, exclude.n = T){

    data.out <- get.boot.dataset(d, site1 = a.site, site2 = b.site, site.1.ref = a.ref, site.2.ref = b.ref, measure)
    boot.out <- get.boot.t2.p.val(data.org = data.out$data.org, data.out$group.ix, B = B, one.ref)
    if(!exclude.n){
        row.out <-  matrix(c(sum(data.out$group.ix == 1), sum(data.out$group.ix == 3), sum(data.out$group.ix == 2), boot.out$statistic, boot.out$pval), nrow = 1)
    }else{
        row.out <-  matrix(c(boot.out$statistic, boot.out$pval), nrow = 1)
    }
    rownames(row.out) <- paste(measure, ', ', a.site, sep = '')
    return(row.out)
}






# this will be used in the one-sample t-tests.
get.ref.means <- function(d, dp, measure = 'bc', ref.ix){
  dp.ref <- dp[ref.ix, ]
  n.ref  <- dim(dp.ref)[1]
  out    <- matrix(NA, n.ref, n.ref)

  if(measure == 'bc'){
    for(i in 1:n.ref){
      for(j in 1:n.ref){
        out[i,j] <- get.bc.mg(dp.ref[i,], dp.ref[j,])
      }
    }
  }else if(measure == 'corr.s'){
    for(i in 1:n.ref){
      for(j in 1:n.ref){
        out[i,j] <- get.corr(dp.ref[i,], dp.ref[j,], method = 'spearman')
      }
    }
  }else if(measure == 'corr.p'){
    for(i in 1:n.ref){
      for(j in 1:n.ref){
        out[i,j] <- get.corr(dp.ref[i,], dp.ref[j,], method = 'pearson')
      }
    }
  }

  return(mean(out[diag(dim(out)[1]) == 0])) # mean excluding the diagonal
}





# this will be used in the one-sample t-tests.
get.ref.stats <- function(d, dp, measure = 'bc', ref.ix, stat = 'mean'){
    dp.ref <- dp[ref.ix, ]
    n.ref  <- dim(dp.ref)[1]
    out    <- matrix(NA, n.ref, n.ref)

    if(measure == 'bc'){
        for(i in 1:n.ref){
            for(j in 1:n.ref){
                out[i,j] <- get.bc.mg(dp.ref[i,], dp.ref[j,])
            }
        }
    }else if(measure == 'corr.s'){
        for(i in 1:n.ref){
            for(j in 1:n.ref){
                out[i,j] <- get.corr(dp.ref[i,], dp.ref[j,], method = 'spearman')
            }
        }
    }else if(measure == 'corr.p'){
        for(i in 1:n.ref){
            for(j in 1:n.ref){
                out[i,j] <- get.corr(dp.ref[i,], dp.ref[j,], method = 'pearson')
            }
        }
    }

    if(stat == 'mean'){
        return(mean(out[diag(dim(out)[1]) == 0])) # mean excluding the diagonal
    }else if(stat == 'median'){
        return(median(out[diag(dim(out)[1]) == 0])) # mean excluding the diagonal
    }else if(stat == 'sd'){
        return(sd(out[diag(dim(out)[1]) == 0])) # mean excluding the diagonal
    }
}






print.nice <- function(x, y, digits = 2){
    out <- paste(round(x, digits), ' (', round(y, digits), ')', sep = '')
    return(out)
}



get.ref.mn.sd.nice <- function(d, dp, measure = 'bc', ref.ix, digits = 2){

    mn <- get.ref.stats(d, dp, measure, ref.ix, stat = 'mean')
    sd <- get.ref.stats(d, dp, measure, ref.ix, stat = 'sd')

    out <- print.nice(mn, sd, digits)
    return(out)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             BOOTSTRAP FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function get.t2
# Expects:
# data     = x1, y1, x2, y2: x = ref 1, y = ref 2, 1 = site 1, 2 site 2,
# group.ix = 1,2,3: 1 = samples from site 1 only, 2 = site 2 only, 3 = site 1 and 2
# boot.ix  = bootstrap index, sampling within groups with replacement

get.boot.t2 <- function(data.null.mx, group.ix, boot.ix = NULL, one.ref = F){
    if(one.ref){
        out <- get.boot.t2.one.ref(data.null.mx, group.ix, boot.ix)
    }else{
        out <- get.boot.t2.two.ref(data.null.mx, group.ix, boot.ix)
    }
}



get.boot.t2.one.ref <- function(data.null.mx, group.ix, boot.ix = NULL){

    if(!is.null(boot.ix)){
        d.mx <- data.null.mx[boot.ix, ]
    }else{
        d.mx <- data.null.mx
    }
    # the number of subjects in each group should be the same for group.ix and group.ix[boot.ix]
    n1 <- sum(group.ix == 1)
    n3 <- sum(group.ix == 3)
    n2 <- sum(group.ix == 2)

    # covariance matrix

    a11 <- 0
    a12 <- 0
    a22 <- 0

    if(n1 > 0){
        a11 <- var(d.mx[,'x1'], na.rm = T)
    }

    if(n2 > 0){
        a22 <- var(d.mx[,'x2'], na.rm = T)
    }

    if(n3 > 0){
        a12 <- cov(x = d.mx[,'x1'], y = d.mx[,'x2'], use = 'pairwise.complete.obs')
    }

    A <- matrix(0, nrow = 2, ncol = 2)
    A[1, 1] <- a11/(n1 + n3)
    A[2, 2] <- a22/(n3 + n2)
    A[1, 2] <- a12*(n3/((n1+n3)*(n3+n2)))
    A[2, 1] <- A[1, 2]

    lambda <- matrix(c(1,-1), byrow = F, ncol = 1)

    D <- matrix(c(mean(d.mx[,'x1'], na.rm = T) - mean(d.mx[,'x2'], na.rm = T)), ncol = 1)

    sigma <- t(lambda)%*%A%*%lambda

    t2 <- t(D)%*%solve(sigma)%*%D

    return(t2)

}


get.boot.t2.two.ref <- function(data.null.mx, group.ix, boot.ix = NULL){

    if(!is.null(boot.ix)){
        d.mx <- data.null.mx[boot.ix, ]
    }else{
        d.mx <- data.null.mx
    }
    # the number of subjects in each group should be the same for group.ix and group.ix[boot.ix]
    n1 <- sum(group.ix == 1)
    n3 <- sum(group.ix == 3)
    n2 <- sum(group.ix == 2)

    # covariance matrix
    a11 <- matrix(0, 2, 2)
    a12 <- matrix(0, 2, 2)
    a22 <- matrix(0, 2, 2)

    if(n1 > 0){
        a11[1,1] <- var(d.mx[,'x1'], na.rm = T)
        a11[1,2] <- var(d.mx[,'x1'], d.mx[,'y1'], na.rm = T)
        a11[2,1] <- a11[1,2]
        a11[2,2] <- var(d.mx[,'y1'], na.rm = T)
    }

    if(n2 > 0){
        a22[1,1] <- var(d.mx[,'x2'], na.rm = T)
        a22[1,2] <- var(d.mx[,'x2'], d.mx[,'y2'], na.rm = T)
        a22[2,1] <- a22[1,2]
        a22[2,2] <- var(d.mx[,'y2'], na.rm = T)
    }

    if(n3 > 0){
        a12[1,1] <- cov(x = d.mx[,'x1'], y = d.mx[,'x2'], use = 'pairwise.complete.obs')
        a12[1,2] <- cov(x = d.mx[,'x1'], y = d.mx[,'y2'], use = 'pairwise.complete.obs')
        a12[2,1] <- cov(x = d.mx[,'y1'], y = d.mx[,'x2'], use = 'pairwise.complete.obs')
        a12[2,2] <- cov(x = d.mx[,'y1'], y = d.mx[,'y2'], use = 'pairwise.complete.obs')
    }

    A <- matrix(0, nrow = 4, ncol = 4)
    A[1:2, 1:2] <- a11/(n1 + n3)
    A[3:4, 3:4] <- a22/(n3 + n2)
    A[1:2, 3:4] <- a12*(n3/((n1+n3)*(n3+n2)))
    A[3:4, 1:2] <- t(a12)*(n3/((n1+n3)*(n3+n2)))

    lambda <- matrix(c(1,0,-1,0,0,1,0,-1), byrow = F, ncol = 2)

    D <- matrix(c(mean(d.mx[,'x1'], na.rm = T) - mean(d.mx[,'x2'], na.rm = T),
                  mean(d.mx[,'y1'], na.rm = T) - mean(d.mx[,'y2'], na.rm = T)),
                ncol = 1)

    sigma <- t(lambda)%*%A%*%lambda

    t2 <- t(D)%*%solve(sigma)%*%D

    return(t2)
}



get.boot.ix <- function(group.ix){
    ix.start <- 1:length(group.ix)
    ix.1 <- group.ix == 1
    ix.3 <- group.ix == 3
    ix.2 <- group.ix == 2

    boot.ix <- sample(ix.start[ix.1], size = sum(ix.1), replace = T)
    boot.ix <- c(boot.ix, sample(ix.start[ix.3], size = sum(ix.3), replace = T))
    boot.ix <- c(boot.ix, sample(ix.start[ix.2], size = sum(ix.2), replace = T))

    return(matrix(boot.ix, ncol = 1))
}



# data = matrix x1, y1, x2, y2: x = ref 1, y = ref 2, 1 = site 1, 2 site 2
# group.ix = 1,2,3: 1 = samples from site 1 only, 2 = site 2 only, 3 = site 1 and 2
get.boot.t2.p.val <- function(data.org, group.ix, B, one.ref = F){
    if(one.ref){
        out <- get.boot.t2.p.val.one.ref(data.org, group.ix, B)
    }else{
        out <- get.boot.t2.p.val.two.ref(data.org, group.ix, B)
    }
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

    group1.ix <- (id %in% u.ids[u.ids.group == 1]) & site1.ix
    group2.ix <- (id %in% u.ids[u.ids.group == 2]) & site2.ix
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





get.row <- function(d, a.site, a.ref, b.site, b.ref, measure, B = 1000, one.ref = F, boot.t.t2 = F, main.paper = F){
    d.a <- d[d$ref == a.ref & d$site == a.site, measure]
    d.b <- d[d$ref == b.ref & d$site == b.site, measure]
    if(boot.t.t2 == F){ # mean, median, IQR
        desc.stats <- c(get.mean.median.iqr(d.a), get.mean.median.iqr(d.b))
    }else if(main.paper == T){
        desc.stats <- c(get.mean.sd(d.a), get.mean.sd(d.b))
    }else{ # mean, median
        desc.stats <- c(get.mean.median(d.a), get.mean.median(d.b))
    }
    out <- matrix(c(desc.stats,
                    t.test.wrap(d, a.site, a.ref, b.site, b.ref, measure),
                    t.test.jack(d, a.site, a.ref, b.site, b.ref, measure),
                    get.boot.row(d, a.site, a.ref, b.site, b.ref, measure, B, one.ref, boot.t.t2)), nrow = 1)
    if(main.paper & boot.t.t2 == T){
        out <- matrix(out[c(1:4, 13:15, 5:12)], nrow = 1)
    }
    rownames(out) <- paste(a.site, ' vs ', b.site, sep = '')
    return(out)
}







get.boot.row <- function(d, a.site, a.ref = T, b.site, b.ref = F, measure, B = 1000, one.ref = F, boot.t.t2 = F){
    data.out <- get.boot.dataset(d, site1 = a.site, site2 = b.site, site.1.ref = a.ref, site.2.ref = b.ref, measure)
    if(boot.t.t2 == F){ # just t2
        boot.out <- get.boot.t2.p.val(data.org = data.out$data.org, data.out$group.ix, B = B, one.ref)
        row.out <- c(boot.out$statistic, boot.out$pval, sum(data.out$group.ix == 1), sum(data.out$group.ix == 3), sum(data.out$group.ix == 2))
    }else if(boot.t.t2 == T & one.ref == T){ # t-test, t2
        boot.out <- get.boot.t.test.t2.p.val.one.ref(data.org = data.out$data.org, data.out$group.ix, B = B)
        row.out <- c(boot.out$t.statistic, boot.out$t.pval, boot.out$t2.statistic, boot.out$t2.pval, sum(data.out$group.ix == 1), sum(data.out$group.ix == 3), sum(data.out$group.ix == 2))
    }else{
        print('no such option in get.boot.row')
    }
    return(row.out)
}





# d = full dataset, including distances
# a.site, b.site = text, sites for a and b
# a.ref, b.ref = T/F
# measure = text, either 'bc' or 'corr'
get.boot.t.test.stat <- function(data.null.mx, group.ix, boot.ix = NULL){

    if(!is.null(boot.ix)){
        d.mx <- data.null.mx[boot.ix, ]
    }else{
        d.mx <- data.null.mx
    }
    t.stat <- t.test(na.omit(d.mx[,1]), na.omit(d.mx[,2]))$statistic # Welch's t-test statistic
    return(abs(t.stat))
}



# data = matrix x1,  x2: x = ref 1, 1 = site 1, 2 site 2
# group.ix = 1,2,3: 1 = samples from site 1 only, 2 = site 2 only, 3 = site 1 and 2
get.boot.t.test.t2.p.val.one.ref <- function(data.org, group.ix, B){
    x1.mn <- mean(data.org[,1], na.rm = T) # mean(x1, na.rm = T)
    x2.mn <- mean(data.org[,2], na.rm = T) # mean(x2, na.rm = T)

    x.mn <- mean(c(data.org[,1], data.org[,2]), na.rm = T) # mean(c(x1, x2), na.rm = T)

    mean.mx <- matrix(rep(c(x1.mn, x2.mn), dim(data.org)[1]), nrow =  dim(data.org)[1], ncol = 2, byrow = T)

    data.null.mx <- data.org - mean.mx

    colnames(data.null.mx) <- c('x1', 'x2')

    t.obs.org     <- get.boot.t.test.stat(data.org, group.ix = group.ix, boot.ix = NULL)
    t.obs.org.mx  <- matrix(t.obs.org, nrow = B, ncol = 1)

    t2.obs.org    <- get.boot.t2(data.org, group.ix = group.ix, boot.ix = NULL, one.ref = T)
    t2.obs.org.mx <- matrix(t2.obs.org, nrow = B, ncol = 1)

    t.boot  <- matrix(NA, nrow = B, ncol = 1)
    t2.boot <- matrix(NA, nrow = B, ncol = 1)

    for(i in 1:B){
        boot.ix    <- get.boot.ix(group.ix)
        t.boot[i]  <- get.boot.t.test.stat(data = data.null.mx, group.ix = group.ix, boot.ix = boot.ix)
        t2.boot[i] <- get.boot.t2(data = data.null.mx, group.ix = group.ix, boot.ix = boot.ix, one.ref = T)
    }

    t.pval <- mean(t.boot >= t.obs.org.mx)
    t2.pval <- mean(t2.boot >= t2.obs.org.mx) # wrt original data
    return(list(t.statistic = t.obs.org, t.pval = t.pval, t2.statistic = t2.obs.org, t2.pval = t2.pval))
}



