find.ruleSubsets = function(recips, combs.df) { # Removes rule combinations which are subsets of another combination
  if (nrow(recips) == 0) { return(rep(F, nrow(combs.df))) }
  
  dup.rule = apply(recips, 1, function(r) {
    r1 = apply(combs.df, 1, function(x) all(combs.df[r[1],] %in% x, na.rm=T))        
    r2 = apply(combs.df, 1, function(x) all(combs.df[r[2],] %in% x, na.rm=T))
    r1&r2
  })
  
  matrixStats::rowAnys(dup.rule)
}

find.recipM = function(M, combs, ppm, ppm.mass) { # Removes rule combinations that are reciporical.  Eg. +H -H
  
  ms = rowSums(array(M$m[combs], dim = dim(combs)), na.rm=T)
  zs = rowSums(array(M$z[combs], dim = dim(combs)), na.rm=T)
  
  zs.e = outer(zs, zs, "+")
  ms.e = outer(ms, ms, "+")
  
  zs.e[upper.tri(zs.e, T)] = NA
  ms.e[upper.tri(zs.e, T)] = NA
  
  which(zs.e == 0 & abs(ms.e) / ppm.mass * 1E6 < ppm, arr.ind=T)
}

build.mCombs = function(M, n, ppm, ppm.mass) {
 combs = lapply(seq(n), gtools::combinations, n = length(M$m), repeats.allowed = T)

 combs.df = do.call(plyr::rbind.fill.matrix, combs)
 
 recips = find.recipM(M, combs.df, ppm, ppm.mass)
 
 badcombs = find.ruleSubsets(recips, combs.df)
 
 combs.df[!badcombs,,drop=F]
}

