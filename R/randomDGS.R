randomDGS <- function (n.DGS = 1, exo.facs = c(""), seed.1 = NULL, seed.2 = NULL) {
 
  if (is.null(seed.1) | is.null(seed.2)) {
   
      errmsg <- paste0("For reproducability, please specify 'seed.1' and 'seed.2'.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  exl <- length(exo.facs)
  
  set.seed(seed.1)
  outcome.columns <- matrix(sample(c(0,1), size = n.DGS*2^exl, replace = TRUE), 
                            ncol = n.DGS)
  
  if (any(apply(outcome.columns, 2, function (x) mean(x) %in% c(0,1)))) {
   
      errmsg <- paste0("At least one model is TRUE (all output function values
                        are 1) or FALSE (all output function values are 0). Please 
                        respecify 'seed.1'.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  # create the configuration matrix for each truth table and add the outcome Z
  dat <- expand.grid(rep(list(c(0,1)), exl))
  dat <- cbind(dat, outcome.columns)
  colnames(dat) <- c(exo.facs, paste0("OUT", seq(n.DGS)))
 
  # set up list that includes all models for each solution
  mods <- vector("list", n.DGS)
 
  for (i in seq(n.DGS)) {
  
       mods[[i]] <- eQMC(dat, exo.facs = exo.facs, outcome = colnames(dat[exl + i]))$solution
  }
 
  nomps <- sapply(mods, length) # number of models per solution
 
  # if there are structural ambiguities, randomly sample one model from the
  # full model space of the solution
 
  set.seed(seed.2)
  for (i in seq(n.DGS)) {
 
       mods[[i]] <- mods[[i]][sample.int(nomps[i], 1)]
  }
 
  mods.dat <- list(DGS = unlist(lapply(unlist(mods, recursive = FALSE), paste0, collapse = " + ")), 
                   tt = dat)
 
  return(mods.dat)
}