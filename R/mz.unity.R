#' Perform a combinatorial search for peak relationships in mass spectral data.
#' 
#' \code{mz.unity.search} searches <m,z> values supplied in \code{B} and \code{M} for combinations which sum to a <m,z> value in \code{A}.
#' 
#' This function encapsulates the combinatorial search for simple and complex peak relationships.  The search is very general and can be tailored to a specific type of relationship with the parameters \code{BM.limits} and \code{M}.  \code{M} allows you to specify additional granular formula which are not present in A and B yet could relate peaks.  This inlcudes things like isotopic masses and small ions.  \code{BM.limits} specifies the combinatorial depth with which to search.  For isotopes in which no intermediate masses are missing BM.limits could be \code{1, 1, 1}.  We can search beyond one missing peak with \code{1,2,1}.
#' 
#' Several sets of granular formula are included by default.  They include isotopes (\code{M.iso}), charge carriers (\code{M.z}), and neutral formula (\code{M.n}).  These can be combined or extended to build a search for any set of relationships. It is generally useful to keep the set of granular formula in small groups appropriate for the combinatorial depth.  For example, isotopes generally require only a single step, whereas we may want to search for cross-polarity relationships which require the exchange of two or more charge carriers.
#' 
#' @param A Matrix containing columns named "m" and "z".  Column "m" is an exact mass. Column "z" is an integer charge.
#' @param B Matrix containing columns named "m", "z", and "d".   Column "d" is one of {-1, 0, 1} denoting weather this species may be added (1) or subtracted (-1) or both (0).
#' @param M Matrix. See B.B and M differ in that they can be included in different proportions as specified in variable \code{BM.limits}.
#' @param ppm Numeric. The mass error allowed when defining mass equality.  This is per included species in A and B.
#' @param BM.limits Matrix containing columns names "M.min", "M.max", "B.n".  This defines the combinatorial depth with which to seach B and M.  For the number of species specified in B.n we search the range of values specified in \code{M.min:M.max}
#' 
#' @return A matrix. Rows correspond to a detected relationship.  Columns include "ppm", "A", multiple "B.*", and multiple "M.*".  A, B, and M are integers referencing the rows of the supplied matrices.  ppm is the mass error when comparing the sum of B and M with the mass of A.
#' 
#' @seealso See \url{https://github.com/nathaniel-mahieu/mz-unity} for examples. \code{\link[mz-unity.graph]} \code{\link[mz-unity.specgraph]}
#' 
#' @export
#' 
mz.unity.search = function(A, B, M, ppm, BM.limits = cbind(M.min = c(0, 0, 0), M.max = c(0, 1, 3), B.n = c(1, 2, 3))) {
  
  M2 = do.call(rbind, list(
    M[which(M[,"d"] == 0),c("m", "z")],
    M[which(M[,"d"] == 0),c("m", "z")] * -1,
    M[which(M[,"d"] == -1),c("m", "z")] * -1,
    M[which(M[,"d"] == 1),c("m", "z")]
  ))
  Mids = c(
    which(M[,"d"] == 0),
    which(M[,"d"] == 0) * -1,
    which(M[,"d"] == -1) * -1,
    which(M[,"d"] == 1)
  )
  
  B2 = do.call(rbind, list(
    B[which(B[,"d"] == 0),c("m", "z")],
    B[which(B[,"d"] == 0),c("m", "z")] * -1,
    B[which(B[,"d"] == -1),c("m", "z")] * -1,
    B[which(B[,"d"] == 1),c("m", "z")]
  ))
  Bids = c(
    which(B[,"d"] == 0),
    which(B[,"d"] == 0) * -1,
    which(B[,"d"] == -1) * -1,
    which(B[,"d"] == 1)
  )
  
  
  rels = execute.search(as.data.frame(A), as.data.frame(B2), as.data.frame(M2), ppm, BM.limits)
  
  rels = as.matrix(rels)
  
  rels[,grep('^M$|^M.', colnames(rels))] = Mids[rels[,grep('^M$|^M.', colnames(rels))]]
  rels[,grep('^B$|^B.', colnames(rels))] = Bids[rels[,grep('^B$|^B.', colnames(rels))]]
  
  as.data.frame(rels)
  }


execute.search = function(A, B, M, ppm, BM.limits) {
  matches = list()

  for (i in seq_along(BM.limits[,"B.n"])) {
    x = BM.limits[i,"B.n"]
    
    B.combs = gtools::combinations(length(B$m), x, repeats.allowed=T)
    Bm.sums = rowSums(array(B$m[B.combs], dim = dim(B.combs)))
    Bz.sums = rowSums(array(B$z[B.combs], dim = dim(B.combs)))
 
    ABm.d = outer(A$m, Bm.sums, "-")
    ABm.s = outer(A$m, Bm.sums, "+")
    ABz.d = outer(A$z, Bz.sums, "-")

    mats = M.search(ABm.d, ABm.s, ABz.d, M, A, ppm, M.limits = BM.limits[i,])
    
    mats = lapply(mats, function(l) {
      ABm.d.ind = R.utils::arrayIndex(l$d, dim(ABm.d))
      
      l$A = ABm.d.ind[,1]
      l$B = B.combs[ABm.d.ind[,2], ,drop=F]
      l
      })
  
    matches = c(matches, mats)
  }
  

  templist = suppressWarnings(lapply(matches, do.call, what=data.frame))
  templist = lapply(templist, function(x) {
    c = colnames(x)
    c[c == "M"] = "M.1"
    c[c == "B"] = "B.1"
    
    colnames(x) = c
    x[,c!="d"]
    })
  
  df = do.call(plyr::rbind.fill, templist)
  
  colns = grep('^M$|^M.', names(df))
  df[,colns][as.matrix(df[,colns]) == 0] = NA
  
  colns = grep('^M$|^M.|^B$|^B.', names(df))
  selfs = rowSums(!is.na(df[,colns])) <= 1
  
  df = df[!selfs,,drop=F]
  
  
  bs = grep('^B$|^B.', colnames(df))
  dfconcerns = which(rowSums(!is.na(df[,bs,drop=F])) == 1)
  
  if (length(dfconcerns) < 2) {return(df)}
  
  df2 = df[dfconcerns, ,drop=F]
  
  revs = sapply(seq(nrow(df2)), function(i) {
    x = df2[i,]

    if (sum(!is.na(x[bs])) > 1) { return(F) }
    
    if (any( which(df2[i:nrow(df2),"A"] %in% x["B.1"]) %in% which(df2[i:nrow(df2),"B.1"] %in% x["A"]))) {
      return (T)
    }
    return(F)
    
    })
  
  df = df[!(seq(nrow(df)) %in% dfconcerns[revs]), ,drop=F]

  df[,c(
    grep('^A$', colnames(df)),
    grep('^ppm$', colnames(df)),
    grep('^M', colnames(df)),
    grep('^B', colnames(df))
    ),drop=F]
  }



M.search = function(d, s, zd, M, A, ppm, M.limits) {

  d.which = which(
    matrixStats::rowAnys(outer(c(d), c(outer(max(M$m), seq(from = M.limits["M.min"], to = M.limits["M.max"]), "*"))+1, "<")) & 
      matrixStats::rowAnys(outer(c(d), c(outer(min(M$m), seq(from = M.limits["M.min"], to = M.limits["M.max"]), "*"))-1, ">"))
    )
  
  M.combs.all = build.mCombs(M, M.limits["M.max"], ppm, 100)
  dimnames(M.combs.all) = NULL
  M.combs.n = rowSums(!is.na(M.combs.all))

  mat.l = list()  
  for (y in M.limits["M.min"]:M.limits["M.max"]) {

    M.combs = M.combs.all[M.combs.n == y, ,drop=F]
    sM   = rowSums(array(M$m[M.combs], dim = dim(M.combs)), na.rm=T)
    sM.z = rowSums(array(M$z[M.combs], dim = dim(M.combs)), na.rm=T)
    if (all(is.na(sM))) { M.combs = array(0, dim=c(1,1)); sM = 0; sM.z = 0 }
    
    ppms = outer(d[d.which], sM, "-") / s[d.which] * 1E6
    mats = which(abs(ppms) < ppm & outer(zd[d.which], sM.z, "-") == 0, arr.ind=T)
    
    if (is.null(nrow(mats))) {mats = matrix(numeric(), nrow=0, ncol=2)}
    
    mat.l = c(
      mat.l,
      list(list(
        d = d.which[mats[,1]],
        M = M.combs[mats[,2],,drop=F],
        ppm = ppms[mats]
        )
      ))
  }
  
  mat.l
}



#' Transverse a graph collecting all connected nodes.
#' 
#' \code{aggregate.self} returns a vector of group assignments containing nodes whose parents already belong to the group.  This is convenient for gathering compound groups.
#' 
#' @param rels Matrix. Representation of the relationship graph with column names "A" and "B.*".
#'
#' @return Numeric. A matrix containing group assignments for the correspond columns in "A" and "B.*"
#' 
#' @export aggregate.self
#' 
aggregate.self = function(rels, force = F) {
  
  pMer = rels
  pMer = as.matrix(pMer)
  
  As = pMer[,grep('^A', colnames(pMer))]
  Bs = pMer[,grep('^B', colnames(pMer)), drop=F]
  
  groups = list()
  names = unique(c(As, Bs))
  
  repeat {
    n = names[1]
    
    repeat {
      prev.len = length(n)
      
      n = unique(c(
        n, 
        As[
          matrixStats::rowAlls(Bs %in% c(n, NA), dim = dim(Bs))
          ], 
        Bs[
          As %in% n & (rowSums(!is.na(Bs)) == 1) & !is.na(Bs)
          ]
      ))
      
      if (length(n) == prev.len) break();
    }
    
    groups = c(groups, list(n))
    names = names[-which(names %in% n)]
    if (length(names) == 0) break();
  }
  
  # Discards subsets of a larger group
  subgroups = sapply(groups, function(x) {
    sum(sapply(groups, function(y) { all(x %in% y) })) > 1
  })
  groups = groups[!subgroups]
  
  
  mer.groups = groups
  group.is = unlist(sapply(seq(groups), function(x) {rep(x, length(groups[[x]])) }))
  names(group.is) = unlist(groups)
  
  # Check that each peak belongs to one group
  ps = unique(unlist(groups))
  ns = sapply(ps, function(x) {
    which(sapply(groups, function(y) {x %in% y}))
  })
  test.conflict = which(sapply(ns, length) > 1)
  if (length(test.conflict) > 0 & !force) { 
    stop(paste(length(test.conflict), "Mers belong to conflicting groups.", paste(test.conflict, collapse=","))) 
  } else if (length(test.conflict) > 0 & force) {
      cons = ps[test.conflict]
      group.is[as.character(cons)] = 0
  }
  
  

  
  abcols = pMer[,grep('^A|^B', colnames(pMer))]
  group.mat = array(group.is[as.character(abcols)], dim=dim(abcols))
  colnames(group.mat) = colnames(abcols)
  
  group.mat
}


#' Transverse a directed graph collecting all child nodes.
#' 
#' \code{aggregate.children} returns a vector of nodes which are children at all depths of a specified node.  This is convenient for gathering isotopic spectra.
#' 
#' @param rels Matrix. Representation of the relationship graph with column names "A" and "B.*".
#' @param node Integer. The node at which to start looking for children.
#'
#' @return Numeric. A vector containing all child nodes.
#' 
#' @export aggregate.children
#' 
aggregate.children = function(rels, node) {
  
  B = node
  
  pMer = rels
  pMer = as.matrix(pMer)
  
  As = pMer[,grep('^A', colnames(pMer))]
  Bs = pMer[,grep('^B', colnames(pMer)), drop=F]
  
  names = unique(c(As, Bs))
  
  n = B
  repeat {
    prev.len = length(n)
    
    n = unique(c(
      n, 
      As[
        matrixStats::rowAlls(Bs %in% c(n, NA), dim = dim(Bs))
        ]
    ))
    
    if (length(n) == prev.len) break();
  }
  
  n
}
