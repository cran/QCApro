implicantMatrix <- function (noflevels, raw = FALSE, arrange = FALSE) {
  
 aEmat <- mintermMatrix(noflevels + 1)
  
  if (arrange) {
      aEmat <- sortMatrix(aEmat)
      sum.zeros <- apply(aEmat, 1, function(idx) sum(idx == 0))
      aEmat <- aEmat[order(sum.zeros, decreasing=TRUE), ]
  }
  
 return(structure(list(aE = aEmat - 1, raw = raw), class = "aE"))

}

