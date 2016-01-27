eQMCLoop <- function(data, outcome = "", neg.out = FALSE, exo.facs = c(""), n.cut = 1,
                     incl.cut1 = 1, incl.cut0 = 1, minimize = c("1"), sol.type = "ps", 
                     row.dom = FALSE, min.dis = FALSE, omit = c(), dir.exp = c(), 
                     details = FALSE, show.cases = FALSE, use.tilde = FALSE, 
                     use.letters = FALSE, inf.test = c(""),
                     relation = "suf", ...) {
    
    check.object <- verify.mqca(data, outcome, exo.facs)
    
    exo.facs <- check.object$exo.facs
    
    data <- data[, unique(exo.facs, check.object$outcome)]
    
    eQMC.list <- solution.list <- vector(mode = "list", length = length(outcome))
    names(eQMC.list) <- names(solution.list) <- outcome
    
    for (i in seq(length(outcome))) {
        
     exo.facs <- names(data)[-which(names(data) == check.object$outcome[i])]
        
        eQMC.list[[i]] <- eQMC(data, outcome = outcome[i], exo.facs = exo.facs, n.cut = n.cut, 
                               incl.cut1 = incl.cut1, incl.cut0 = incl.cut0, neg.out = neg.out, 
                               minimize = minimize, sol.type = sol.type, row.dom = row.dom, 
                               min.dis = min.dis, relation = relation, 
                               show.cases = show.cases, 
                               ... = ...)
    }
    
    return(structure(eQMC.list, class = "mqca"))
    
}
