eQMC <- function(data, outcome = c(""), neg.out = FALSE, exo.facs = c(""), 
                 relation = "suf", n.cut = 1, incl.cut1 = 1, incl.cut0 = 1, 
                 minimize = c("1"), sol.type = "ps", row.dom = FALSE, min.dis = FALSE, 
                 omit = c(), dir.exp = c(), details = FALSE, show.cases = FALSE, 
                 inf.test = c(""), use.tilde = FALSE, use.letters = FALSE, ...) {
    
  # specify solution types (argument 'include' to be replaced in future version)
  if (sol.type == "cs") {include <- c("")}
  else if (sol.type == "cs+") {include <- c("C")}
  else if (sol.type == "ps") {include <- c("?")}
  else if (sol.type == "ps+") {include <- c("?", "C")}
  else {errmsg <- paste0("The specified solution type, ", sol.type,", does not 
                          exist. Either \"ps\" or \"ps+\" should normally be used.")
        cat("\n")
        stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
    
  metacall <- match.call()
  other.args <- list(...)
  PRI <- FALSE

  if ("PRI" %in% names(other.args)) {
  
      if (is.logical(other.args$PRI)) {
   
          PRI <- other.args$PRI
      }
  }
    
  # before: print.truth.table <- details & !is.tt(data)
  print.truth.table <- FALSE
  
  if (all(include == "")) {
        
      if (!is.null(dir.exp)) {
        
       errmsg <- paste0("Directional expectations cannot be specified in 
                         conjunction with the conservative solution. Please use
                         sol.type = \"ps\" or sol.type = \"ps+\" for this purpose.")
       cat("\n")
       stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
        
      else {
        
          include <- minimize
      }
  }
    
  if (!is.tt(data)) {
     
      # if "matrix", first turn into "data.frame" 
      if (is.matrix(data)) {
      
          data <- as.data.frame(data)
      }
        
      if (length(outcome) > 1) {
            
          return(eQMCLoop(data = data, outcome = outcome, neg.out = neg.out, 
                          exo.facs = exo.facs, n.cut = n.cut, incl.cut1 = incl.cut1, 
                          incl.cut0 = incl.cut0, minimize = minimize, include = include, 
                          row.dom = row.dom, min.dis = min.dis, omit = omit, 
                          dir.exp = dir.exp, details = details, show.cases = show.cases,
                          use.tilde = use.tilde, use.letters = use.letters, 
                          inf.test = inf.test, relation = relation, 
                          ... = ...))
      }
   
      verify.qca(minimize = minimize)
      
      outcome.copy <- outcome
      indata <- data # important before altering the outcome, if multi-value
      
      if (grepl("[{]", outcome)) {
       
          outcome <- substr(outcome, 1, regexpr("[{]", outcome)[1] - 1)
          outcome.value <- as.numeric(substr(outcome.copy, regexpr("[{]", outcome.copy)[1] + 1, 
                                                           regexpr("}", outcome.copy)[1] - 1))
       
          data[, outcome] <- ifelse(data[, outcome] == outcome.value, 1, 0)
      }
      
      complete <- FALSE
      if ("complete" %in% names(other.args)) {
       
       complete <- other.args$complete
      }
      
      tt <- truthTable(data = data, outcome = outcome, neg.out = neg.out,
                       exo.facs = exo.facs, n.cut = n.cut, incl.cut1 = incl.cut1, 
                       incl.cut0 = incl.cut0, complete = complete, 
                       show.cases = show.cases, use.letters = use.letters,
                       PRI = PRI)
      
      # if no exogenous factors are specified, use all factors except that of the outcome
      if (all(gsub("\\s", "", exo.facs, perl = TRUE) == "")) {
       
          exo.facs <- colnames(data)[-which(colnames(data) == outcome)]
      }
      
      # make sure all relevant labels are uppercase
      idx.names <- match(colnames(data[, c(exo.facs, outcome)]), colnames(data))
      colnames(data)[idx.names] <- toupper(colnames(data)[idx.names])
      exo.facs <- toupper(exo.facs)
      outcome <- toupper(outcome)
        
      data <- data[, c(exo.facs, outcome)]
        
      tt$initial.data <- indata
      indata <- data # data is already altered in outcome value, if initially multi-value
      recdata <- tt$recoded.data
      colnames(indata) <- c(exo.facs, outcome)
      
      # dir.exp should now be a list, in the same order as exo.facs
      dir.exp <- verify.dir.exp(recdata, outcome, exo.facs, dir.exp)
        
      if (!is.null(dir.exp)) {
        
          names(dir.exp) <- toupper(names(dir.exp))
      }
        
      rowsNotMissing <- which(tt$tt$OUT != "?")
  }
    
  else { # data already is a tt
    
      chminimize <- c(0, 1)[which(0:1 %in% minimize)]
      chinclude <- c(0, 1)[which(0:1 %in% include)]
        
      if (length(chinclude) > 0) {
        
          if (any(chinclude != chminimize)) {
        
              chinclude <- chinclude[which(chinclude != chminimize)]
              stop(paste("\nYou cannot include ", chinclude, " since you want to minimize ", chminimize, ".\n\n", sep=""), call. = FALSE)
          }
      }
    
      # check if minimize has both 1 and 0
        if (length(chminimize) == 2) {
        
            stop("\nPositive and negative minterms must not be combined.\n\n", call. = FALSE)
        }
        
        tt <- data
        indata <- tt$initial.data
        recdata <- tt$recoded.data
        exo.facs <- colnames(recdata)[seq(length(tt$noflevels))]
        outcome <- colnames(recdata)[ncol(recdata)]
        use.letters <- tt$opts$use.letters
        
        rowsNotMissing <- which(tt$tt$OUT != "?")
        if (any(tt$tt$OUT == "?")) {
            missings <- which(tt$tt$OUT == "?")
            tt$tt <- tt$tt[-missings, ]
        }
        
        neg.out <- tt$neg.out
        
        dir.exp <- verify.dir.exp(recdata, outcome, exo.facs, dir.exp)
        if (!is.null(dir.exp)) {
            names(dir.exp) <- toupper(names(dir.exp))
        }
    }
    
    
    uplow <- TRUE
    noflevels <- tt$noflevels
     # check if the column names are not already letters
    alreadyletters <- sum(nchar(colnames(recdata)[-ncol(recdata)])) == ncol(recdata) - 1
    
    output <- list()
    output$tt <- tt
    output$opts$print.truth.table <- print.truth.table
    
    tt$tt[, seq(length(exo.facs))] <- as.data.frame(lapply(tt$tt[, seq(length(exo.facs))], function(x) {
        x[x %in% c("-", "dc")] <- -1
        return(as.numeric(x))
    }))
    
    expl.incl <- unique(c(minimize, include)) # here "include" may contain contradictions; missings are irrelevant as they were already erased
    subset.tt <- tt$tt[, "OUT"] %in% expl.incl
    expl.matrix <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    expl.matrix <- matrix(as.numeric(expl.matrix), ncol = length(noflevels)) + 1
    rownames(expl.matrix) <- tt$indexes[subset.tt]
    
    subset.tt <- !tt$tt[, "OUT"] %in% expl.incl
    excl.matrix <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    excl.matrix <- matrix(as.numeric(excl.matrix), ncol = length(noflevels)) + 1
    
    subset.tt <- tt$tt[, "OUT"] %in% minimize
    
    if (all(!subset.tt)) {
        stop(paste0("\nNone of the output function values is minimized. Please check the truth table.\n\n"), call. = FALSE)
    }
    
    inputt <- as.matrix(tt$tt[subset.tt, seq(length(noflevels))])
    rownms <- rownames(inputt)
    inputt <- matrix(as.numeric(inputt), ncol=length(noflevels)) + 1
    inputcases <- tt$cases[rowsNotMissing][subset.tt]
    
    nofcases1 <- sum(tt$tt$n[tt$tt$OUT == 1])
    nofcases0 <- sum(tt$tt$n[tt$tt$OUT == 0])
    nofcasesC <- sum(tt$tt$n[tt$tt$OUT == "C"])
    
    tomit <- logical(nrow(expl.matrix))
    tomitinputt <- logical(nrow(inputt))
    if (is.matrix(omit)) {
        cnoflevels <- noflevels
        for (i in seq(ncol(omit))) {
            if (any(omit[, i] < 0)) {
                omit[, i][omit[, i] < 0] <- noflevels[i]
                cnoflevels[i] <- noflevels[i] + 1
            }
        }
        omitrows <- drop(rev(c(1, cumprod(rev(cnoflevels))))[-1] %*% t(omit)) + 1
        tomit <- rownames(expl.matrix) %in% omitrows
        tomitinputt <- rownms %in% omitrows
        excl.matrix <- rbind(excl.matrix, omit + 1)
    }
    else if (is.vector(omit)) {
        tomit <- rownames(expl.matrix) %in% omit
        tomitinputt <- rownms %in% omit
        excl.matrix <- unique(rbind(excl.matrix, getRow(noflevels, as.numeric(omit)) + 1))
    }
    
    output$excluded <- sort(drop(rev(c(1, cumprod(rev(noflevels))))[-1] %*% t(excl.matrix - 1)) + 1)
    expl.matrix <- expl.matrix[!tomit, , drop = FALSE]
    inputt <- inputt[!tomitinputt, , drop = FALSE]
    inputcases <- inputcases[!tomitinputt]
    
    if (nrow(expl.matrix) == 0) {
        
        stop("\nAll observed minterms have the same function value. Please check the truth table.\n\n", call. = FALSE)
    }
    
    incl.rem <- "?" %in% include
    
    if (nrow(excl.matrix) == 0 & incl.rem) {
    
        stop("\nAll observed minterms have the same function value. Please check the truth table.\n\n", call. = FALSE)
    }
    
     # expl.matrix needs to be unaltered for the incl.rem argument
    expressions <- expl.matrix
    
    recdata[, exo.facs] <- as.data.frame(lapply(recdata[, exo.facs], function(x) {
        x[x %in% c("-", "?", "dc")] <- -1
        return(as.numeric(x))
    }))
    
     # check if the data has multiple values
    if (any(recdata[, seq(ncol(recdata) - 1)] > 1)) {
        uplow <- FALSE
        use.tilde <- FALSE
    }
    
    if (use.tilde) {
        uplow <- FALSE
    }
    
    collapse <- ifelse(alreadyletters & uplow | use.tilde, "", "*")
    changed <- FALSE
    
    
     # if not already letters and user specifies using letters for conditions, change it
    if (use.letters & !alreadyletters) {
        colnames(expressions) <- colnames(inputt) <- colnames(expl.matrix) <- LETTERS[seq(ncol(inputt))]
        changed <- TRUE
        collapse <- ifelse(!uplow | use.tilde, "*", "")
    }
    else {
        colnames(expressions) <- colnames(inputt) <- colnames(expl.matrix) <- colnames(recdata[, seq(ncol(inputt))])
        if (use.tilde) {
            collapse <- "*"
        }
    }
    
    output$initials <- writePrimeimp(inputt, collapse=collapse, uplow=uplow, use.tilde=use.tilde)
    initial <- drop(rev(c(1, cumprod(rev(noflevels))))[-1] %*% t(inputt - 1)) + 1
    
    minExpressions <- function(expressions) {
        minimized <- TRUE
        while (any(minimized)) {
            minimized <- logical(nrow(expressions))
            distance <- dist(expressions, method="manhattan")
            distance <- as.matrix(distance)
            distance[!upper.tri(distance)] <- NA
            tbc <- as.matrix(which(distance == 1, arr.ind=TRUE))
            
            if (nrow(tbc) > 0) {
                differences <- t(apply(tbc, 1, function(idx) expressions[idx[1], ] != expressions[idx[2], ]))
                result <- matrix(0, nrow=0, ncol=ncol(differences))
                for (i in seq.int(nrow(differences))) {
                    stable.values <- expressions[tbc[i, 1], !differences[i, , drop=FALSE], drop=FALSE]
                    subset.minimize <- apply(expressions[, !differences[i, , drop=FALSE], drop=FALSE], 1, function(x) all(x == stable.values))
                    if (sum(subset.minimize) == noflevels[differences[i, ]]) {
                        minimized[subset.minimize] <- TRUE
                        minimization.result <- expressions[tbc[i, 1], , drop=FALSE]
                        minimization.result[differences[i, ]] <- 0
                        result <- rbind(result, as.vector(minimization.result))
                    }
                }
            }
            
            if (sum(minimized) > 0) {
                expressions <- rbind(expressions[!minimized, ], unique(result))
            }
        }
        return(expressions)
    }
    
    expressions <- minExpressions(expressions)
    
    c.sol <- p.sol <- getSolution(expressions = expressions, collapse = collapse,
                                  uplow = uplow, use.tilde = use.tilde, inputt = inputt,
                                  row.dom = row.dom, initial = initial, min.dis = min.dis,
                                  ... = ...)
    
    mbase <- rev(c(1, cumprod(rev(noflevels + 1))))[-1]
    
    if (incl.rem) {
        
            expressions <- sort(setdiff(findSupersets(noflevels + 1, expl.matrix), findSupersets(noflevels + 1, excl.matrix)))
            expressions <- .Call("removeRedundants", expressions, noflevels, mbase, PACKAGE="QCApro")
        
        expressions <- getRow(noflevels + 1, expressions)
        
        colnames(expressions) <- colnames(inputt)
        
        p.sol <- getSolution(expressions = expressions, collapse = collapse, 
                             uplow = uplow, use.letters = use.letters, use.tilde = use.tilde, inputt = inputt, 
                             row.dom = row.dom, initial = initial, min.dis = min.dis, 
                             ... = ...)
        
    }
    
    output$PIs <- p.sol$all.PIs
    output$PIchart <- structure(list(p.sol$mtrx), class = "pic")
    output$primes <- p.sol$reduced$expressions
    output$solution <- p.sol$solution.list[[1]]
    output$essential <- p.sol$solution.list[[2]]
    
        expr.cases <- rep(NA, nrow(p.sol$reduced$expressions))
            
            tt.rows <- createString(inputt - 1, collapse=collapse, uplow, use.tilde)
            
            mtrxlines <- demoChart(rownames(p.sol$reduced$expressions), tt.rows, ifelse((use.letters & uplow) | (alreadyletters & uplow), "", "*"))
            
            for (l in seq(length(expr.cases))) {
                expr.cases[l] <- paste(inputcases[mtrxlines[l, ]], collapse="; ")
            }
        
        if (length(p.sol$solution.list[[1]]) == 1) {
            
            listIC <- pof(p.sol$reduced$expressions - 1, tt$outcome, indata, showc=TRUE, cases=expr.cases, neg.out=neg.out,
                          relation = "sufficiency", via.eQMC = TRUE)
            listIC$opts$show.cases <- show.cases
        }
        else {
            
            listIC <- pof(p.sol$reduced$expressions - 1, tt$outcome, indata, showc=TRUE, cases=expr.cases, neg.out=neg.out,
                          relation = "sufficiency", via.eQMC = TRUE, solution.list=output$solution, essential=output$essential)
            listIC$opts$show.cases <- show.cases
        }
        
        output$pims <- listIC$pims
        
        listIC$pims <- NULL
        output$IC <- listIC
    
    
    output$numbers <- c(OUT1 = nofcases1, OUT0 = nofcases0, OUTC = nofcasesC, 
                        Total = nofcases1 + nofcases0 + nofcasesC)
    
    output$inputcases <- inputcases
    
    output$opts$minimize <- minimize
    output$opts$neg.out <- neg.out
    output$opts$details <- details
    output$opts$show.cases <- show.cases
    output$opts$use.letters <- use.letters
    output$opts$collapse <- collapse
    
    if (PRI) {
        output$opts$PRI <- other.args$PRI[1]
    }
    
    output$SA <- lapply(p.sol$solution.list[[1]], function(x) {
        p.expressions <- p.sol$reduced$expressions[x, , drop=FALSE]
        
        temp <- apply(p.expressions, 1, function(pr) {
            indices <- rev(which(!pr))
            
            SA <- NULL
            for (k in indices) {
                if (is.null(SA)) {
                    SA <- drop(mbase %*% pr) + sum(mbase[!pr])
                }
                tempSA <- SA
                for (lev in seq(noflevels[k] - 1)) {
                    tempSA <- c(tempSA, SA + mbase[k]*lev)
                }
                SA <- tempSA
            }
            return(SA)
        })
        
        if (all(is.null(temp))) {
            return(NULL)
        }
        else {
            temp <- sort(unique(as.vector(unlist(temp))))
            temp <- temp[!temp %in% drop(mbase %*% t(inputt))]
            if (length(temp) > 0) {
                SA <- getRow(noflevels + 1,  temp + 1) - 1
                colnames(SA) <- colnames(inputt)
                rownames(SA) <- drop(c(rev(cumprod(rev(noflevels))), 1)[-1] %*% t(SA)) + 1
                return(SA)
            }
            else {
                return(NULL)
            }
        }
    })
    
    
    prettyNums <- formatC(seq(length(p.sol$solution.list[[1]])), digits = nchar(length(p.sol$solution.list[[1]])) - 1, flag = 0)
    names(output$SA) <- paste0("M", prettyNums)
    
    if (!is.null(dir.exp) & all(include != c(""))) {
        
        i.sol <- vector("list", length(c.sol$solution.list[[1]])*length(p.sol$solution.list[[1]]))
        index <- 1
        
        for (c.s in seq(length(c.sol$solution.list[[1]]))) {
            
            c.expressions <- c.sol$reduced$expressions[c.sol$solution.list[[1]][[c.s]], , drop = FALSE]
            
            for (p.s in seq(length(p.sol$solution.list[[1]]))) {
                
                p.expressions <- p.sol$reduced$expressions[p.sol$solution.list[[1]][[p.s]], , drop = FALSE]
                
                dir.exp.matrix <- matrix(matrix(ncol = length(exo.facs), nrow = 0))
                
                for (i in seq(nrow(c.expressions))) {
                    comp <- c.expressions[i, ]
                    
                    for (j in seq(nrow(p.expressions))) {
                        pars <- p.expressions[j, ]
                        
                        dir.exp.temp <- rep(-1, length(pars))
                        equals <- comp[pars > 0] == pars[pars > 0]
                        
                        if (all(equals) > 0) {
                            res <- lapply(dir.exp, function(x) return(-1))
                            equals <- which(pars > 0)
                            for (k in equals) {
                                res[[k]] <- sort(drop(as.numeric(pars[k] - 1)))
                            }
                            
                            dir.exp.temp[equals] <- c.expressions[i, equals] - 1
                            notequals <- setdiff(which(comp > 0), equals)
                            if (length(notequals) > 0) {
                                for (k in notequals) {
                                    
                                    dir.exp.k <- unique(c(names(
                                      dir.exp[[exo.facs[k]]]
                                    )[dir.exp[[exo.facs[k]]] == 1], c.expressions[i, k] - 1))
                                    
                                    if (length(dir.exp.k) != noflevels[k]) {
                                        equals <- sort(c(equals, k))
                                        res[[k]] <- sort(drop(as.numeric(dir.exp.k)))
                                    }
                                }
                            }
                            
                            dir.exp.matrix <- rbind(dir.exp.matrix, expand.grid(res))
                            
                        }
                        else {
                        }
                    }
                }
                
                names(i.sol)[index] <- paste0("C", c.s, "P", p.s)
                
                EC <- matrix(ncol = length(exo.facs), nrow = 0)
                colnames(EC) <- colnames(inputt)
                
                for (dir.exp.i in seq(nrow(dir.exp.matrix))) {
                    dir.exp.x <- dir.exp.matrix[dir.exp.i, ]
                    
                    if (!is.null(output$SA[[p.s]])) {
                        SArows <- apply(output$SA[[p.s]], 1, function(x) {
                            return(all(x[dir.exp.x >= 0] == dir.exp.x[dir.exp.x >= 0]))
                        })
                        
                        subSA <- output$SA[[p.s]][SArows, , drop = FALSE]
                        EC <- rbind(EC, subSA[setdiff(rownames(subSA), rownames(EC)), , drop = FALSE])
                    } 
                }
                
                i.sol[[index]]$EC <- EC[order(as.numeric(rownames(EC))), , drop = FALSE]
                i.sol[[index]]$DC <- output$SA[[p.s]][setdiff(rownames(output$SA[[p.s]]),
                  rownames(EC)), , drop = FALSE]
                i.sol[[index]]$NSEC <- matrix(ncol = ncol(EC), nrow = 0)
                colnames(i.sol[[index]]$NSEC) <- colnames(EC)
                
                nsecs <- TRUE
                
                while (nsecs) {
                    expl.matrix.i.sol <- unique(rbind(expl.matrix, i.sol[[index]]$EC + 1))
                    
                    tomit <- logical(nrow(expl.matrix.i.sol))
                    if (is.matrix(omit)) {
                        cnoflevels <- noflevels
                        for (i in seq(ncol(omit))) {
                            if (any(omit[, i] < 0)) {
                                omit[, i][omit[, i] < 0] <- noflevels[i]
                                cnoflevels[i] <- noflevels[i] + 1
                            }
                        }
                        tomit <- rownames(expl.matrix) %in% (drop(rev(c(1, cumprod(rev(cnoflevels))))[-1] %*% t(omit)) + 1)
                    }
                    else if (is.vector(omit)) {
                        if (is.numeric(omit)) {
                            tomit <- rownames(expl.matrix) %in% as.character(omit)
                        }
                    }
                    
                    expl.matrix.i.sol <- expl.matrix.i.sol[!tomit, , drop = FALSE]
                    
                    expressions <- minExpressions(expl.matrix.i.sol)
                    
                    i.sol.index <- getSolution(expressions = expressions,
                      collapse = collapse, uplow = uplow, use.tilde = use.tilde,
                      inputt = inputt, row.dom = row.dom, initial = initial,
                      min.dis = min.dis, ... = ...)
                    
                    i.sol.index$expressions <- i.sol.index$expressions[rowSums(i.sol.index$mtrx) > 0, , drop = FALSE]
                    
                    if (nrow(i.sol[[index]]$EC) > 0) {
                        nsec <- !vector(length = nrow(i.sol[[index]]$EC))
                        
                        for (i in seq(nrow(i.sol.index$expressions))) {
                            i.sol.PI <- i.sol.index$expressions[i, ]
                            for (j in seq(length(nsec))) {
                                j.EC <- i.sol[[index]]$EC[j, ]
                                
                                if (all(i.sol.PI[i.sol.PI > 0] == (j.EC[i.sol.PI > 0] + 1))) {
                                    nsec[j] <- FALSE
                                }
                            }
                        }
                        
                        nsecs <- any(nsec)
                    }
                    else {
                        nsecs <- FALSE
                    }
                    
                    if (nsecs) {
                        i.sol[[index]]$NSEC <- rbind(i.sol[[index]]$NSEC,
                          i.sol[[index]]$EC[which(nsec), , drop = FALSE])
                        i.sol[[index]]$EC <- i.sol[[index]]$EC[-which(nsec), , drop = FALSE]
                    }
                }
                
                i.sol[[index]]$PIchart <- structure(list(i.sol.index$mtrx), class = "pic")
                i.sol[[index]]$c.sol <- c.sol$solution.list[[1]][[c.s]]
                i.sol[[index]]$p.sol <- p.sol$solution.list[[1]][[p.s]]
                i.sol[[index]]$solution <- i.sol.index$solution.list[[1]]
                i.sol[[index]]$essential <- i.sol.index$solution.list[[2]]
                i.sol[[index]]$primes <- i.sol.index$reduced$expressions
                
                expr.cases <- rep(NA, nrow(i.sol.index$reduced$expressions))
                    
                    tt.rows <- createString(inputt - 1, collapse=collapse, uplow, use.tilde)
                    
                    mtrxlines <- demoChart(rownames(i.sol.index$reduced$expressions), tt.rows, ifelse((use.letters & uplow) | (alreadyletters & uplow), "", "*"))
                    
                    for (l in seq(length(expr.cases))) {
                        expr.cases[l] <- paste(inputcases[mtrxlines[l, ]], collapse="; ")
                    }
                
                if (length(i.sol.index$solution.list[[1]]) == 1) {
                    i.sol[[index]]$IC <- pof(i.sol.index$reduced$expressions - 1,
                      outcome, indata, showc = TRUE, cases = expr.cases,
                      relation = "sufficiency", neg.out = neg.out,
                      via.eQMC = TRUE)
                    
                    i.sol[[index]]$IC$opts$show.cases <- show.cases
                }
                else {
                    i.sol[[index]]$IC <- pof(i.sol.index$reduced$expressions - 1,
                      outcome, indata, showc = TRUE, cases = expr.cases,
                      relation = "sufficiency", neg.out = neg.out,
                      via.eQMC = TRUE, solution.list = i.sol.index$solution.list[[1]],
                      essential = i.sol.index$solution.list[[2]])
                    i.sol[[index]]$IC$opts$show.cases <- show.cases
                }
                i.sol[[index]]$pims <- i.sol[[index]]$IC$pims
                i.sol[[index]]$IC$pims <- NULL
                index <- index + 1
            }
        }
        output$i.sol <- i.sol
    }
    
    # put SAs and ECs in data frames
    
    output$SA <- lapply(output$SA, as.data.frame)
    
    if (any(names(output) == "i.sol")) {
        for (i in seq(length(output$i.sol))) {
            output$i.sol[[i]]$EC <- as.data.frame(output$i.sol[[i]]$EC)
        }
    }
    
    if (!is.tt(data)) {
        output$tt$outcome <- outcome.copy
    }
    
    output$relation <- relation
    output$call <- metacall
    
    return(structure(output, class = "qca"))
}
