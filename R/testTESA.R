testTESA <- function (data, outcome = "", neg.out = FALSE, exo.facs = c(""), 
                      n.cut = 1, incl.cut1 = 1, incl.cut0 = 1) {
  
  # if no exogenous factors are specified, use all factors except that of the outcome
  if (all(gsub("\\s", "", exo.facs, perl = TRUE) == "")) {
      
      exo.facs <- colnames(data)[-which(colnames(data) == outcome)]
  }
 
  # create the full truth table
  tt.b <- truthTable(data, outcome = outcome, neg.out = neg.out, exo.facs = exo.facs, 
                     n.cut = 1, incl.cut1 = incl.cut1, incl.cut0 = incl.cut0, 
                     complete = TRUE)
  
  tt <- tt.b$tt[tt.b$tt$OUT == "?", seq(match("OUT", colnames(tt.b$tt)))]
  
  names.tt <- names(tt)
  
  for (i in seq(length(names.tt))) {
  
       tt[, i] <- ifelse(tt[, i] == 1, names.tt[i], tolower(names.tt[i]))
  }
 
  # identify all minimally necessary conditions
  necs <- superSubset(data, outcome = outcome, neg.out = neg.out, 
                      exo.facs = exo.facs, incl.cut = incl.cut1)
  
  necs <- as.list(strsplit(names(necs$coms), "[+]"))
  
  necs.test <- apply(do.call("rbind", lapply(necs, function (x) {
    apply(tt, 1, function (y) all(x %in% y))})), 2, any)
  
  SA.list <- lapply(eQMC(tt.b)$SA, row.names)
  
  model.test <- sapply(SA.list, function (x) {
  
    round(sum(x %in% names(which(necs.test == TRUE)))/length(x), 3) * 100}
  )
 
  return(model.test)
}