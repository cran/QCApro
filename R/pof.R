pof <- function(setms, outcome, data, neg.out = FALSE, relation = "suf",
         inf.test = c(""), incl.cut1 = 0.75, incl.cut0 = 0.5, ...) {
    
  funargs <- as.list(match.call())
    
  other.args <- list(...)
    
  recursive <- "recursive" %in% names(other.args)
  via.eQMC <- "via.eQMC" %in% names(other.args)
  force.rows <- "force.rows" %in% names(other.args)
    
  if (recursive) {
      mins <- other.args$mins
      outcome <- other.args$vo
      sum.outcome <- other.args$so
      pims <- other.args$pims
      incl.cov <- matrix(NA, nrow = ncol(mins), ncol = 4)
  }
    
  else {
      
      outcomename <- ""
        
      if (all(is.character(outcome)) & length(outcome) == 1) {
      
          if (missing(data)) {
              cat("\n")
              stop("The data argument is missing, with no default.\n\n", call. = FALSE)
          }
          
          else {
             
              if (grepl("[{]", outcome)) {
                 
                  outcome.copy <- outcome
                  outcome <- substr(outcome, 1, regexpr("[{]", outcome)[1] - 1)
                  outcome.value <- as.numeric(substr(outcome.copy, regexpr("[{]", outcome.copy)[1] + 1, 
                                                     regexpr("}", outcome.copy)[1] - 1))
                  
                  data[, outcome] <- ifelse(data[, outcome] == outcome.value, 1, 0)
              }
                
              outcomename <- toupper(outcome)
              outcome <- data[, toupper(outcome)]
          }
      }
      
      else if (is.vector(outcome)) {
      
          outcomename <- as.character(funargs$outcome)
            
          if (outcomename[1] == "[") {
          
              condnumber <- suppressWarnings(as.numeric(outcomename[length(outcomename)]))
        
              if (is.na(condnumber)) {
                 
                  outcomename <- toupper(outcomename[length(outcomename)])
              }
              
              else {
              
                  outcomename <- toupper(colnames(get(outcomename[2]))[condnumber])
              }
          }
            
          else if (outcomename[1] == "$") {
                
              outcomename <- toupper(outcomename[length(outcomename)])
          }
      }
      
      else {
          
          cat("\n")
          stop("The outcome should be either a column name in the data or a vector of values.\n\n", call. = FALSE)
      }
        
      if (!(relation %in% c("necessity", "sufficiency", "suf", "nec"))) {
          
          cat("\n")
          stop("The relationship should be either \"necessity\" or \"sufficiency\".\n\n", call. = FALSE)
      }
        
      if (!missing(data)) {
            
          if (length(outcome) != nrow(data)) {
        
              cat("\n")
              stop("The outcome's length should be the same as the number of rows in the data.\n\n", call. = FALSE)
          }
            
          if (any(outcomename %in% names(data))) {
           
              noflevels <- truthTable(data, outcome=outcomename, via.pof=TRUE)
          }
          
          else {
          
              cat("\n")
              stop("The outcome was not found in the data.\n\n", call. = FALSE)
          }
            
          colnames(data) <- toupper(colnames(data))
          outcomename <- toupper(outcomename)
          exo.facs <- colnames(data)[-which(colnames(data) == outcomename)]
          data <- data[, c(exo.facs, outcomename)]
      }
        
        pims <- FALSE
        
        if (is.data.frame(setms)) {
            if (missing(outcome)) {
                cat("\n")
                stop("Outcome is missing, with no default.\n\n", call. = FALSE)
            }
            else {
                # colnames(setms) <- toupper(colnames(setms))
                exo.facs <- colnames(setms)
                
                if (missing(data)) { # outcome was already checked to be (or coerced to) a vector
                    if (nrow(setms) == length(outcome)) {
                        newdata <- cbind(setms, outcome)
                        colnames(newdata)[ncol(newdata)] <- outcomename
                        pims <- TRUE
                    }
                    else {
                        cat("\n")
                        stop("The length of outcome should be the same as the number of rows in \"setms\".\n\n", call. = FALSE)
                    }
                }
                else {
                    if (nrow(setms) == nrow(data)) {
                        data <- cbind(setms, outcome)
                        colnames(data)[ncol(data)] <- outcomename
                        pims <- TRUE
                    }
                    else {
                        cat("\n")
                        stop("The number of rows in \"setms\" should be the same as the number of rows in the data.\n\n", call. = FALSE)
                    }
                }
            }
        }
        else if (is.matrix(setms)) {
            if (missing(data)) {
                cat("\n")
                stop("The data argument is missing, with no default.\n\n", call. = FALSE)
            }
            
            if (ncol(setms) == length(exo.facs)) {
                setms[setms < 0] <- -1
                setms <- setms + 1
            }
            else {
                cat("\n")
                stop("The number of columns in the \"setms\" does not match the number of conditions.\n\n", call. = FALSE)
            }
        }
        else if (is.vector(setms)) {
            setms <- suppressWarnings(as.numeric(setms))
            setms <- setms[!is.na(setms)]
            
            if (length(setms) == 0) {
                cat("\n")
                stop("The \"setms\" argument does not contain any numbers.\n\n", call. = FALSE)
            }
            
            if (force.rows) {
                if (missing(data)) {
                    cat("\n")
                    stop("The data argument is missing, with no default.\n\n", call. = FALSE)
                }
                
                if (any(table(setms) > 1) | any(setms == 0)) {
                    cat("\n")
                    stop("The \"setms\" argument does not appear to be a vector of row numbers.\n\n", call. = FALSE)
                }
                
                setms <- getRow(noflevels + 1, setms)
                
            }
            else {
                if (length(setms) == length(outcome)) {
                    if (all(setms >= 0 & setms <= 1) | any(table(setms) > 1)) {
                        newdata <- cbind(setms, outcome)
                        exo.facs <- as.character(funargs$setms)
                        
                        if (length(exo.facs) > 1) {
                            if (exo.facs[1] == "[" || exo.facs[1] == "[[") {
                                condnumber <- suppressWarnings(as.numeric(exo.facs[length(exo.facs)]))
                                
                                if (is.na(condnumber)) {
                                    exo.facs <- exo.facs[length(exo.facs)]
                                }
                                else {
                                    obj <- exo.facs[2]
                                    objdollar <- unlist(strsplit(obj, split="\\$"))
                                    if (length(objdollar) > 1) {
                                        obj <- get(objdollar[1])
                                        for (i in seq(2, length(objdollar))) {
                                            obj <- obj[[objdollar[i]]]
                                        }
                                    }
                                    else {
                                        obj <- get(obj)
                                    }
                                    exo.facs <- colnames(obj)[condnumber]
                                }
                            }
                            else if (exo.facs[1] == "$") {
                                exo.facs <- exo.facs[length(exo.facs)]
                            }
                            else if (exo.facs[1] == "-") {
                                exo.facs <- exo.facs[length(exo.facs)]
                                condollar <- unlist(strsplit(exo.facs, split="\\$"))
                                
                                
                                if (length(condollar) > 1) {
                                    obj <- get(condollar[1])
                                    lastcond <- condollar[length(condollar)]
                                    
                                    condsplit <- unlist(strsplit(lastcond, split=""))
                                    
                                    if (condsplit[length(condsplit)] == "]") {
                                        condindex <- unlist(strsplit(lastcond, split="\\["))
                                        
                                        condollar[length(condollar)] <- condindex[1]
                                        
                                        for (i in seq(2, length(condollar))) {
                                            obj <- obj[[condollar[i]]]
                                        }
                                        
                                        condindex <- condindex[length(condindex)]
                                        
                                        condindex <- gsub(",", "", condindex)
                                        condindex <- gsub("\\]", "", condindex)
                                        condindex <- gsub("^[[:space:]]+|[[:space:]]+$", "", condindex)
                                        condindex <- gsub("\"", "", condindex)
                                        
                                        condnumber <- suppressWarnings(as.numeric(condindex))
                                        
                                        if (is.na(condnumber)) {
                                            exo.facs <- condindex
                                        }
                                        else {
                                            exo.facs <- colnames(obj)[condnumber]
                                        }
                                    }
                                    else {
                                        exo.facs <- lastcond
                                    }
                                }
                                else {
                                    condsplit <- unlist(strsplit(condollar, split=""))
                                    if (condsplit[length(condsplit)] == "]") {
                                        condindex <- unlist(strsplit(condollar, split="\\["))
                                        obj <- get(condindex[1])
                                        condindex <- condindex[length(condindex)]
                                        
                                        condindex <- gsub(",", "", condindex)
                                        condindex <- gsub("\\]", "", condindex)
                                        condindex <- gsub("^[[:space:]]+|[[:space:]]+$", "", condindex)
                                        condindex <- gsub("\"", "", condindex)
                                        
                                        condnumber <- suppressWarnings(as.numeric(condindex))
                                        
                                        if (is.na(condnumber)) {
                                            exo.facs <- condindex
                                        }
                                        else {
                                            exo.facs <- colnames(obj)[condnumber]
                                        }
                                    }
                                }
                            }
                            else {
                                cat("\n")
                                stop("Unknown type of input in the \"setms\" argument.\n\n", call. = FALSE)
                            }
                        }
                        
                        colnames(newdata) <- c(exo.facs, outcomename)
                        pims <- TRUE
                    }
                    else {
                        if (any(setms > 5)) {
                            cat("\n")
                            stop("Assuming this is a vector of row numbers, data argument is missing with no default (try force.rows = TRUE).\n\n", call. = FALSE)
                        }
                    }
                }
                else {
                    if (missing(data)) {
                        cat("\n")
                        stop("Data argument is missing, or the length of \"setms\" is not equal to the length of outcome.\n\n", call. = FALSE)
                    }
                    else {
                        setms <- getRow(noflevels + 1, setms)
                    }
                }
            }
        }
        else {
            cat("\n")
            stop("The \"setms\" argument is not correct.\n\n", call. = FALSE)
        }
        
        if (missing(data)) {
            data <- as.data.frame(newdata)
            noflevels <- truthTable(data, outcome=outcomename, via.pof=TRUE)
        }
        
        if (is.matrix(setms)) { # necessary here and not above because setms might be a vector and then transformed into a matrix via getRow()
            if (is.null(colnames(setms))) {
                colnames(setms) <- exo.facs
            }
            
            if (is.null(rownames(setms))) {
                use.tilde <- FALSE
                if ("use.tilde" %in% names(other.args)) {
                    rownames(setms) <- writePrimeimp(setms, uplow=all(noflevels == 2), use.tilde=other.args$use.tilde)
                }
                else {
                    rownames(setms) <- writePrimeimp(setms, uplow=all(noflevels == 2))
                }
            }
        }
        
        
        hastime <- logical(length(exo.facs))
        for (i in seq(length(exo.facs))) {
            if (any(data[, i] %in% c("-", "dc", "?"))) {
                hastime[i] <- TRUE
            }
        }
        
        
        if (!pims) {
            setms <- setms[, !hastime, drop=FALSE]
        }
        
        data[, which(hastime)] <- NULL
        exo.facs <- exo.facs[!hastime]
        
        if (neg.out) {
            outcome <- 1 - outcome
        }
        sum.outcome <- sum(outcome)
        
        if (pims) {
            mins <- setms
            
            if (is.vector(setms)) {
                length.expr <- 1
            }
            else {
                length.expr <- ncol(mins)
            }
            incl.cov <- matrix(NA, nrow=length.expr, ncol=4)
        }
        else {
            
            fc <- apply(data[, exo.facs], 2, function(x) any(x %% 1 > 0))
            incl.cov <- matrix(NA, nrow=nrow(setms), ncol=4)
            
            length.expr <- nrow(setms)
            
            mins <- apply(setms, 1, function(e) {
                apply(data[, exo.facs, drop=FALSE], 1, function(v) {
                    
                    if (any(ox <- e[fc] == 1)) {
                        v[fc][ox] <- 1 - v[fc][ox]
                    }
                    
                    if (length(cp <- v[!fc]) > 0) {
                        v[!fc][e[!fc] != cp + 1] <- 0
                        v[!fc][e[!fc] == cp + 1] <- 1
                    }
                    
                    return(min(v[e != 0]))
                })
            })
            
            if (!is.matrix(mins)) { ## bug fix 10.03.2014, if the data contains a single combination, mins is not a matrix but a vector
                mins <- t(as.matrix(mins))
                rownames(mins) <- rownames(data)
            }
        }
    }
    
    
    if (is.vector(mins)) {
        mins <- as.data.frame(mins)
        colnames(mins) <- exo.facs
    }
    
    
    rownames(incl.cov) <- colnames(mins)
    colnames(incl.cov) <- c("incl", "PRI", "cov.r", "cov.u")
    
    pmins <- apply(mins, 2, pmin, outcome)
    primins <- apply(mins, 2, function(x) pmin(x, 1 - outcome, outcome))
    
    if (relation %in% c("necessity", "nec")) {
        primins <- apply(mins, 2, function(x) pmin(x, 1 - x, outcome))
    }
    
    if (!is.matrix(pmins)) { ## bug fix 10.03.2014, if the data contains a single combination, pmins is not a matrix but a vector
        pmins <- t(as.matrix(pmins))
        rownames(pmins) <- rownames(mins)
                ## probably the very same thing happens to primins
        primins <- t(as.matrix(primins))
        rownames(primins) <- rownames(mins)
    }
    
    
    incl.cov[, 1] <- colSums(pmins)/colSums(mins)
    incl.cov[, 2] <- (colSums(pmins) - colSums(primins))/(colSums(mins) - colSums(primins))
    incl.cov[, 3] <- colSums(pmins)/sum.outcome
    
    
    if (relation %in% c("necessity", "nec")) {
        incl.cov[, 1] <- colSums(pmins)/sum.outcome
        incl.cov[, 2] <- (colSums(pmins) - colSums(primins))/(sum.outcome - colSums(primins))
        incl.cov[, 3] <- colSums(pmins)/colSums(mins)
    }
    
    
    maxmins <- fuzzyor(mins) # union
    inclusions <- fuzzyor(pmins)
    prisol <- pmin(maxmins, 1 - outcome, outcome)
    
    if (relation %in% c("necessity", "nec")) {
        prisol <- pmin(maxmins, 1 - maxmins, outcome)
    }
    
    if (ncol(mins) > 1) {
        for (i in seq(nrow(incl.cov))) {
            incl.cov[i, 4] <- incl.cov[i, 3] - sum(pmin(pmins[, i], fuzzyor(pmins[, -i]), outcome))/sum.outcome
        }
    }
    
    
    # solution incl, pri and cov
    sol.incl <- sum(inclusions)/sum(maxmins)
    sol.pri <- (sum(inclusions) - sum(prisol))/(sum(maxmins) - sum(prisol))
    sum.cov <- sum(inclusions)/sum.outcome
    
    
    result.list <- list(incl.cov=as.data.frame(incl.cov, stringsAsFactors = FALSE), relation=relation)
    
    if (!pims & via.eQMC) {
        result.list$sol.incl.cov <- c(incl=sol.incl, PRI=sol.pri, cov=sum.cov)
        result.list$pims <- as.data.frame(mins)
    }
    
    if ("recursive" %in% names(other.args)) {
        return(result.list)
    }
    
    # showc is not a formal argument, therefore is it initiated as FALSE
    showc <- FALSE
    
    if (all(inf.test != "")) {
        verify.inf.test(inf.test, data)
    }
    
    if (inf.test[1] == "binom") {
        
        if ("incl.cut1" %in% names(other.args)) {
            incl.cut1 <- as.numeric(other.args$incl.cut1)
        }
        
        if ("incl.cut0" %in% names(other.args)) {
            incl.cut0 <- as.numeric(other.args$incl.cut0)
        }
        
        if (length(inf.test) > 1) {
            alpha <- as.numeric(inf.test[2]) # already checked if a number between 0 and 1
        }
        else {
            alpha <- 0.05
        }
        
        incl.cov <- as.data.frame(incl.cov, stringsAsFactors = FALSE)
        
        if (relation  %in% c("necessity", "nec")) {
            nofcases <- rep(sum.outcome, ncol(mins))
        }
        else {
            nofcases <- colSums(mins)
        }
        
        success <- as.vector(round(nofcases * as.numeric(incl.cov[, "incl"])))
        
        incl.cov$pval0 <- incl.cov$pval1 <- 0
        
        for (i in seq(length(success))) {
            incl.cov[i, "pval1"] <- binom.test(success[i], nofcases[i], p = incl.cut1, alternative = "less")$p.value
            incl.cov[i, "pval0"] <- binom.test(success[i], nofcases[i], p = incl.cut0, alternative = "greater")$p.value
        }
        
        result.list$incl.cov <- incl.cov
    }
    
    
    if ("showc" %in% names(other.args)) {
        if (other.args$showc) {
            showc <- other.args$showc
            result.list$incl.cov <- cbind(result.list$incl.cov, cases = other.args$cases, stringsAsFactors=FALSE)
        }
    }
    
    if ("solution.list" %in% names(other.args)) {
        solution.list <- other.args$solution.list
        length.solution <- length(solution.list)
        individual <- vector("list", length=length.solution)
        
        for (i in seq(length.solution)) {
            individual[[i]] <- Recall(relation="sufficiency", recursive=TRUE, via.eQMC=TRUE,
                                      mins=mins[, solution.list[[i]], drop=FALSE],
                                      vo = outcome, so = sum.outcome, pims=pims)
        }
        return(structure(list(overall=result.list, individual=individual, essential=other.args$essential, pims=as.data.frame(mins), relation=relation, opts=list(show.cases = showc)), class="pof"))
    }
    else {
        result.list$opts <- list(show.cases = showc)
        return(structure(result.list, class="pof"))
    }
}

