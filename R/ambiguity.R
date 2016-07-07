ambiguity <- function (data, outcome = c(""), neg.out = c(FALSE), 
                       exo.facs = c(""), tuples = c(), incl.cut1 = c(1), 
                       incl.cut0 = c(1), sol.type = c("ps"), row.dom = c(FALSE),
                       min.dis = c(FALSE)) {
 
  if (all(gsub("\\s", "", exo.facs, perl = TRUE) == "")) {
 
      endo.facs <- gsub("[{].{1,}", "", outcome) 
      exo.facs <- colnames(data)[-match(endo.facs, colnames(data))]
  }
 
  verify.ambig(data = data, outcome = outcome, neg.out = neg.out, 
               exo.facs = exo.facs, tuples = tuples, incl.cut1 = incl.cut1, 
               incl.cut0 = incl.cut0, row.dom = row.dom, min.dis = min.dis)
 
  cbn <- n.models <- vector("list", length(tuples))
  names(cbn) <- paste0(tuples, "-tuples")
  names(n.models) <- paste0(tuples, "-tuples")
  
  l.outcome <- length(outcome)
  l.neg.out <- length(neg.out)
  l.tuples <- length(tuples)
  l.incl.cuts <- length(incl.cut1)
  l.sol.type <- length(sol.type)
  l.constraints <- length(row.dom)
 
  for (i in 1:l.tuples) {
       
       cbn[[i]] <- t(combn(exo.facs, tuples[i]))
  }
  
  colnames <- paste0(
     rep(outcome, each = l.neg.out*l.incl.cuts*l.sol.type*l.constraints), ".",
     rep(seq(l.neg.out), each = l.incl.cuts*l.sol.type*l.constraints),
     rep(seq(l.incl.cuts), each = l.sol.type*l.constraints),
     rep(seq(l.incl.cuts), each = l.constraints),
     seq(l.constraints)
  )
  
  for (i in 1:l.tuples) {
  
       n.models[[i]] <- matrix(
         numeric(l.outcome*l.neg.out*l.incl.cuts*l.sol.type*l.constraints*
                 nrow(cbn[[i]])),
         ncol = l.outcome*l.neg.out*l.incl.cuts*l.sol.type*l.constraints,
         dimnames = list(1:nrow(cbn[[i]]), colnames)
       )
  }
  
  for (i in 1:l.tuples) {
  
    for (j in 1:nrow(cbn[[i]])) {
   
      for (m in 1:l.outcome) {
      
        for (p in 1:l.neg.out) {
     
          for (q in 1:l.incl.cuts) {
      
            for (r in 1:l.sol.type) {
       
              for (s in 1:l.constraints)  {
        
                try({sol <- eQMC(data, outcome = outcome[m], neg.out = neg.out[p],
                                 exo.facs = cbn[[i]][j, ], incl.cut1 = incl.cut1[q],
                                 incl.cut0 = incl.cut0[q], sol.type = sol.type[r], 
                                 row.dom = row.dom[s], min.dis = min.dis[s])
                  
       n.models[[i]][j, (m - 1)*(l.neg.out*l.incl.cuts*l.sol.type*l.constraints) +
                        (p - 1)*(l.incl.cuts*l.sol.type*l.constraints) +
                        (q - 1)*(l.sol.type*l.constraints) +
                        (r - 1)*(l.constraints) + s] <- length(sol$solution)
                }, silent = TRUE)
              }
            }
          }
        }
      }
    }
  }
 
  return(list(exo.facs = exo.facs, outcome = outcome, tuples = cbn, 
              n.models = n.models))
}