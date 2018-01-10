superSubset <- function(data, outcome = "", neg.out = FALSE, exo.facs = c(""), 
                        relation = "nec", incl.cut = 1, cov.cut = 0, use.tilde = FALSE, 
                        use.letters = FALSE, ...) {
    
    memcare <- FALSE # to be updated with a future version
    
    other.args <- list(...)
    
    PRI <- FALSE
    if ("PRI" %in% names(other.args)) {
        if (is.logical(other.args$PRI)) {
            PRI <- other.args$PRI[1] # [1] just to make sure only the first value is taken, should someone by mistake provide a vector
        }
    }
    
    incl.cut <- incl.cut - .Machine$double.eps ^ 0.5
    if (cov.cut > 0) {
        cov.cut <- cov.cut - .Machine$double.eps ^ 0.5
    }
    
    if (all(exo.facs == c(""))) {
        exo.facs <- names(data)[-which(names(data) == outcome)]
    }
    
    if (grepl("[{]", outcome)) { # there is a "{" sign in the outcome's name
        outcome <- unlist(strsplit(outcome, split = ""))
        outcome.value <- as.numeric(outcome[which(outcome == "{") + 1L])
        outcome <- paste(outcome[seq(1, which(outcome == "{") - 1L)], collapse="")
        
        if (!any(unique(data[, outcome]) == outcome.value)) {
            cat("\n")
            stop(paste("The value {", outcome.value, "} does not exist in the outcome.\n\n", sep=""), call. = FALSE)
        }
        data[, outcome] <- ifelse(data[, outcome] == outcome.value, 1, 0)
    }
    
    
    verify.data(data, outcome, exo.facs)
    
    if (!relation %in% c("necessity", "sufficiency", "nec", "suf", "sufnec", "necsuf")) {
        stop("\nThe relationship should be either \"necessity\", \"sufficiency\" or \"necsuf\".\n\n", call. = FALSE)
    }
    
    relationcopy <- relation
    if (relation == "sufnec" | relation == "necsuf") {
        cov.cut <- incl.cut
    }
    
    if (relation == "sufnec") {
        relation <- "sufficiency"
    }
    else if (relation == "necsuf") {
        relation <- "necessity"
    }
    
    colnames(data) <- toupper(colnames(data))
    exo.facs <- replacements <- toupper(exo.facs)
    outcome <- toupper(outcome)
    
    data <- data[, c(exo.facs, outcome)]
    nofexo.facs <- length(exo.facs)
    
    
    if (neg.out) {
        data[, outcome] <- 1 - data[, outcome]
    }
    
    uplow <- !use.tilde
    
    fc <- apply(data[, exo.facs], 2, function(x) any(x %% 1 > 0))
    
    if (mv.data <- any(data[, exo.facs] > 1)) {
        uplow <- use.tilde <- FALSE
    }
    
    alreadyletters <- sum(nchar(exo.facs)) == length(exo.facs)
    
    collapse <- ifelse(alreadyletters & uplow | use.tilde, "", "*")
    
    if (use.letters & !alreadyletters) {
        replacements <- LETTERS[seq(length(exo.facs))]
        names(replacements) <- exo.facs
        colnames(data)[seq(length(exo.facs))] <- exo.facs <- replacements
        collapse <- ifelse(!uplow | use.tilde, "*", "")
    }
    
    noflevels <- apply(data[, exo.facs], 2, max) + 1L
    noflevels[fc] <- 2
    mbase <- c(rev(cumprod(rev(noflevels + 1L))), 1)[-1]
    
    
    if (memcare) {
        CMatrix <- .Call("superSubsetMem", as.matrix(data[, exo.facs]), noflevels, mbase, as.numeric(fc), data[, outcome], relation == "necessity", PACKAGE = "QCApro")
    }
    else {
        nk <- mintermMatrix(noflevels + 1L)
        colnames(nk) <- exo.facs
        nk <- nk[-1, ] # first row is always empty
        
        CMatrix <- .Call("superSubset", as.matrix(data[, exo.facs]), nk, as.numeric(fc), data[, outcome], as.numeric(relation == "necessity"), PACKAGE = "QCApro")
    }
    
    
    # to modify this, attributing colnames copies the object and uses too much memory
    expressions <- colnames(CMatrix) <- seq_len(ncol(CMatrix)) + 1L # plus 1 because the first row of the nk matrix was deleted
    lincl <- ifelse(relation %in% c("necessity", "nec"), 2, 1)
    
    expressions <- expressions[CMatrix[lincl, ] >= incl.cut & CMatrix[3 - lincl, ] >= cov.cut]
    
    prev.result <- FALSE
    lexpressions <- length(expressions)
    
    if (lexpressions > 0) {
        if (relation %in% c("sufficiency", "suf")) {
            expressions <- .Call("removeRedundants", expressions, noflevels, mbase, PACKAGE="QCApro")
        }
        
        result.matrix <- getRow(noflevels + 1L, expressions)
        rownames(result.matrix) <- expressions
        colnames(result.matrix) <- exo.facs
        result.matrix <- sortMatrix(result.matrix)
        sum.zeros <- apply(result.matrix, 1, function(idx) sum(idx == 0))
        result.matrix <- result.matrix[order(sum.zeros, decreasing=TRUE), , drop=FALSE]
        #collapsign <- "*"
        row_names <- writePrimeimp(result.matrix, collapse=collapse, uplow=uplow, use.tilde=use.tilde)
        prev.result <- TRUE
        result <- data.frame(incl  = CMatrix[lincl, rownames(result.matrix)],
             PRI   = CMatrix[5, rownames(result.matrix)],
             cov.r = CMatrix[3 - lincl, rownames(result.matrix)],
             stringsAsFactors=FALSE,
             row.names=row_names)
    }
    
    lexprnec <- 0
    if (relation  %in% c("necessity", "nec")) {
        exprnec <- seq_len(ncol(CMatrix)) + 1L
        
        exprnec <- exprnec[CMatrix[4, ] >= incl.cut & CMatrix[3, ] >= cov.cut]
        
        exprnec <- .Call("removeRedundants", exprnec, noflevels, mbase, PACKAGE="QCApro")
        
        exprnec <- setdiff(exprnec, expressions)
        lexprnec <- length(exprnec)
        
        if (lexprnec + lexpressions == 0) {
            cat("\n")
            stop(paste("\nThere are no combinations with incl.cut = ", round(incl.cut, 3), " and cov.cut = ", round(cov.cut, 3), "\n\n", sep=""), call. = FALSE)
        }
        
        if (lexprnec > 0) {
            result.matrix2 <- getRow(noflevels + 1, exprnec)
            rownames(result.matrix2) <- exprnec
            colnames(result.matrix2) <- exo.facs
            result.matrix2 <- sortMatrix(result.matrix2)
            
            sum.zeros <- apply(result.matrix2, 1, function(idx) sum(idx == 0))
            result.matrix2 <- result.matrix2[order(sum.zeros, decreasing=TRUE), , drop=FALSE]
            row_names2 <- writePrimeimp(result.matrix2, collapse="+", uplow=uplow, use.tilde=use.tilde)
            
            if (prev.result) {
                result <- rbind(result, data.frame(incl  = CMatrix[4, rownames(result.matrix2)],
                    PRI   = CMatrix[6, rownames(result.matrix2)],
                    cov.r = CMatrix[3, rownames(result.matrix2)],
                    stringsAsFactors=FALSE,
                    row.names=row_names2))
                row_names <- c(row_names, row_names2)
                result.matrix <- rbind(result.matrix, result.matrix2)
            }
            else {
                result <- data.frame(incl = CMatrix[4, rownames(result.matrix2)],
                    PRI = CMatrix[6, rownames(result.matrix2)],
                    cov.r = CMatrix[3, rownames(result.matrix2)],
                    stringsAsFactors=FALSE,
                    row.names=row_names2)
                row_names <- row_names2
                result.matrix <- result.matrix2
            }
            
        }
    }
    
    if (lexprnec + lexpressions == 0) { # there is no combination which exceeds incl.cut
        cat("\n")
        stop(paste("\nThere are no combinations with incl.cut = ", round(incl.cut, 3), " and cov.cut = ", round(cov.cut, 3), "\n\n", sep=""), call. = FALSE)
    }
    
    
    mins <- matrix(NA, nrow=nrow(data), ncol=nrow(result.matrix))
    for (i in seq(nrow(result.matrix))) {
        mins[, i] <- apply(data[, exo.facs], 1, function(v) {
            e <- result.matrix[i, , drop=FALSE]
            if (any(ox <- e[fc] == 1)) {
                v[fc][ox] <- 1 - v[fc][ox]
            }
            
            if (length(cp <- v[!fc]) > 0) {
                v[!fc][e[!fc] != cp + 1] <- 0
                v[!fc][e[!fc] == cp + 1] <- 1
            }
            if (rownames(e) %in% expressions) {
                return(min(v[e != 0]))
            }
            else {
                return(max(v[e != 0]))
            }
        })
    }
    
    colnames(mins) <- rownames(result)
    rownames(mins) <- rownames(data)
    mins <- as.data.frame(mins)
    
    if (relationcopy == "sufnec") {
        colnames(result) <- c("inclS", "PRI", "inclN")
    }
    else if (relationcopy == "necsuf") {
        colnames(result) <- c("inclN", "PRI", "inclS")
    }
    
    out.list <- list(incl.cov=result, coms=mins, use.letters=use.letters, letters=replacements)
    if (PRI) {
        out.list$PRI <- PRI
    }
    return(out.list)
}



