`solveChart` <-
function(chart, row.dom = FALSE, min.dis = TRUE, ...) {
    
    if (!is.logical(chart)) {
        cat("\n")
        stop("Use a T/F matrix. See demoChart's output.\n\n", call. = FALSE)
    }
    
    other.args <- list(...)
    
    if ("all.sol" %in% names(other.args)) {
        if (is.logical(other.args$all.sol)) {
            min.dis <- !other.args$all.sol
        }
    }
    
    if (!min.dis) {
        row.dom <- FALSE
    }
    
    row.numbers <- seq(nrow(chart))
    
    if (row.dom) {
        row.numbers <- rowDominance(chart)
        chart <- chart[row.numbers, ]
    }
    
    output <- list()
    
    if (all(dim(chart) > 1)) {
         ## initial solution provided by Gabor Grothendieck
         ## the function lp (from package lpSolve) finds a (guaranteed) minimum solution
         # k will be the minimum number of prime implicants necessary to cover all columns
        k <- ceiling(sum(lp("min", rep(1, nrow(chart)), t(chart), ">=", 1)$solution))
        # cat(paste("k: ", k, "\n"))
        
        forceRAM <- 2
        if ("forceRAM" %in% names(other.args)) {
            if (length(other.args$forceRAM) == 1) {
                if (is.numeric(other.args$forceRAM) & other.args$forceRAM > 0) {
                    forceRAM <- other.args$forceRAM
                }
            }
        }
        
        
         # Stop if the matrix with all possible combinations of k PIs has over 2GB of memory
        if ((mem <- nrow(chart)*choose(nrow(chart), k)*8/1024^3) > forceRAM) {
            errmessage <- paste(paste("Too much memory needed (", round(mem, 1), " GB) to solve the PI chart using combinations of", sep=""),
                                   k, "out of", nrow(chart), "minimised PIs, with the PI chart having", ncol(chart), "columns.\n\n")
            cat("\n")
            stop(paste(strwrap(errmessage, exdent = 7), collapse = "\n", sep=""))
        }
        
        if (!min.dis & k < nrow(chart)) {
            
            # if (2^nrow(chart)*2 > .Machine$integer.max) {
            if (nrow(chart) > 29) { # in order to prevent cases where integer.max is larger than 32-bit
                cat("\n")
                stop(paste(strwrap("The PI chart is too large to identify all models.\n\n", exdent = 7), collapse = "\n", sep=""))
            }
            
            output <- .Call("allSol", k, chart*1, PACKAGE="QCApro")
            
            output[output == 0] <- NA
            
        }
        else {
            
             # create a matrix with all possible combinations of k prime implicants
            combos <- combn(nrow(chart), k)
            
             # sol.matrix will be a subset of the chart matrix with all minimum solutions
            output <- combos[, as.logical(.Call("solveChart", t(combos) - 1, chart*1, PACKAGE="QCApro")[[1]]), drop=FALSE]
        }
    }
    else {
        output <- matrix(seq(nrow(chart)))
        
        if (ncol(chart) == 1) {
            output <- t(output)
        }
    } 
    
    # just in case row dominance was applied
    return(matrix(row.numbers[output], nrow=nrow(output)))
}

