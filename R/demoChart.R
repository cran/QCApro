demoChart <- function(primes = c(""), configs = c(""), and.split = "") {
    
  if (and.split != "") {
     
      and.split <- paste0("\\", and.split)
  }
    
  primes.split <- strsplit(primes, and.split)
  configs.split <- strsplit(configs, and.split)
    
  mtrx <- matrix(FALSE, nrow = length(primes), ncol = length(configs))
    
  for (i in seq(nrow(mtrx))) {
   
       for (j in seq(ncol(mtrx))) {
   
            mtrx[i, j] <- all(primes.split[[i]] %in% configs.split[[j]])
       }
  }
    
  colnames(mtrx) <- configs
  rownames(mtrx) <- primes
  return(mtrx)
}

