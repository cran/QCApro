randomDGS <- function (n.DGS = 1, exo.facs = c(""), seed.1 = NULL, seed.2 = NULL, 
                       prob = 0.5, diversity = 1, delete.trivial = FALSE) {
  
  if (is.null(seed.1) | is.null(seed.2)) {
    
      errmsg <- paste0("For reproducability, please specify 'seed.1' and 'seed.2'.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  if (prob <= 0 | prob >= 1) {
    
      errmsg <- paste0("The value to the argument 'prob' must be larger than 0 
                     and smaller than 1. It is currently set to ", prob, ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  if (diversity <= 0 | diversity > 1) {
    
      errmsg <- paste0("The value to the argument 'diversity' must be larger than 0 
                     and not larger than 1. It is currently set to ", diversity, ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  exl <- length(exo.facs)
  
  set.seed(seed.1)
  outcome.cols <- matrix(sample(c(0,1), size = n.DGS*2^exl, replace = TRUE, 
                                prob = c(1 - prob, prob)), 
                         ncol = n.DGS)
  
  if (any(colMeans(outcome.cols) %in% c(0,1))) {
    
      if (delete.trivial == FALSE) {
      
          errmsg <- paste0("At least one model is TRUE (all output function values
                           are 1) or FALSE (all output function values are 0). Please 
                           re-specify 'seed.1', adjust the values to the arguments 
                           'prob' and/or 'diversity', or increase the number of 
                           exogenous factors.")
          cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
    
      else {
      
          trivial <- apply(outcome.cols, 2, function (x) mean(x) %in% c(0,1))
          outcome.cols <- outcome.cols[ , !trivial, drop = FALSE]  
      
          wrnmsg <- paste0(sum(trivial), " trivial structure(s) (TRUE or FALSE) was/were eliminated.")
          cat("\n")
          warning(paste(strwrap(wrnmsg, exdent = 7), collapse = "\n"), 
                  immediate. = TRUE, call. = FALSE)
    }
  }
  
  if (ncol(outcome.cols) == 0) {
    
    errmsg <- paste0("All structures have been trivial (TRUE or FALSE). Please re-specify
                        seed.1, adjust the values to the arguments 'prob' and/or 'diversity',
                        or increase the number of exogenous factors.")
    cat("\n")
    stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  # create minterm matrix for each truth table and add output function values
  dat <- expand.grid(rep(list(c(0,1)), exl))
  dat <- cbind(dat, outcome.cols)
  n.out.cols <- ncol(outcome.cols)
  colnames(dat) <- c(exo.facs, paste0("OUT", seq(n.out.cols)))
  
  # create limited empirical diversity
  
  if (diversity < 1) {
    
    dat <- dat[sample.int(nrow(dat), as.integer(diversity * nrow(dat))),]
  }
  
  # set up list that includes all models for each solution
  mods <- vector("list", n.DGS)
  
  for (i in seq(n.out.cols)) {
    
    mods[[i]] <- eQMC(dat, exo.facs = exo.facs, outcome = colnames(dat[exl + i]))$solution
  }
  
  nomps <- sapply(mods, length) # number of models per solution
  
  # if there are structural ambiguities, randomly sample one model from the
  # full model space of the solution
  
  set.seed(seed.2)
  for (i in seq(n.out.cols)) {
    
    mods[[i]] <- mods[[i]][sample.int(nomps[i], 1)]
  }
  
  mods.dat <- list(DGS = unlist(lapply(unlist(mods, recursive = FALSE), paste0, collapse = " + ")), 
                   tt = dat)
  
  return(mods.dat)
}