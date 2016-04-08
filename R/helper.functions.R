# Input Checking
checkMsZs = function(mszs) {
  if (any(!c("m", "z") %in% colnames(mszs))) { stop("Supplied mass and charges must have columns with names 'm' and 'z'.") }
  if (class(mszs) != 'data.frame') { stop("Supplied mass and charges must be a matrix.") }
}
checkMzs = function(mzs) {
  if (class(mzs) != 'numeric') { stop("Supplied mass and charges must be a numeric vector.") }
}
checkZHypos = function(z) {
  if (class(mzs) != 'numeric') { stop("Supplied charge hypotheses must be a numeric vector.") }
  if (any(sign(z.hypos) < 1)) { stop("Supplied charge hypotheses (z.hypos) must be positive integers representing number of charges.") }
}


mzHypoZ = function(mzs, zs, ds) {
  if (any(zs < 0 | zs%%1!=0)) { stop("Charge sign is mass specific and set in mzs, zs must be positive integers representing number of starting charges.") }
  
  dim = c(length(mzs), length(zs))
  
  df = as.data.frame(cbind(
    m = c(array(abs(mzs), dim=dim)),
    z = c(outer(sign(mzs), zs)),
    n = c(array(seq_along(mzs), dim=dim))
  ))
  
  df$d = ds
  
  df$m = df$m * abs(df$z)
  
  df
}



t01 = function(x) { x[x == 0] = 1; x }

#' Takes a matrix and reindexes cells reminescent of match()
#' 
#' \code{reindex} takes a matrix and reindexes cells reminescent of match()
#' 
#' This is useful for performing a search with a subset of granular formula, and reindexing to some master ID index.
#' 
#' @param pConv Matrix. With integer columns referenced in col.regex to be reindexed.
#' @param col.regex Character. Regex matching to the column names of columns to be reindexed.
#' @param indices Numeric. Numeric vector with new values.  \code{indices[pConv[1,1]]} will be the replacement value
#' 
#' @return A matrix with referenced cells replaced by new values
#' 
#' @export
#' 
reindex = function(pConv, col.regex, indices) {
  if (nrow(pConv) == 0) {return(pConv)}
  
  colns = grep(col.regex, names(pConv))
  is = unname(unlist(pConv[,colns]))
  pConv[,colns] = sign(is)*indices[abs(is)]
  pConv
}

#' Turns a mz.unity relationship output into a graph structured data.frame.
#' @export
expandGraph = function(mat, cols.from = '^B$|^B.', cols.to = '^A$|^A.') {
  from.cols = grep(cols.from, names(mat))
  to.cols = grep(cols.to, names(mat))
  
  names = c()
  pairs = lapply(seq(nrow(mat)), function(i) {
    a = mat[i,to.cols]; a = a[!is.na(a)]
    b = mat[i,from.cols]; b = b[!is.na(b)]
    
    
    r = expand.grid(a,b)
    names <<- c(names, rep(i, nrow(r)))
    
    r
  })
  
  r = do.call(rbind, pairs)
  r$row = names
  r
}