truthTable <- function(data, outcome = "", neg.out = FALSE, exo.facs = c(""), 
                       n.cut = 1, incl.cut1 = 1, incl.cut0 = 1, complete = FALSE, 
                       show.cases = FALSE, sort.by = c(""), decreasing = TRUE, 
                       inf.test = c(""), use.letters = FALSE, ...) {
    
    memcare <- FALSE # to be updated with a future version
    metacall <- match.call()
    other.args <- list(...)
    via.pof <- "via.pof" %in% names(other.args)
    
    if (memcare) {
        complete <- FALSE
    }
    
    # first turn "matrix" into "data.frame" 
    if (is.matrix(data)) {
        
     data <- as.data.frame(data)
    }
    
    # call the verification function for truth tables
    if (!via.pof) {
     
        verify.tt(data = data, outcome = outcome, exo.facs = exo.facs, 
                  complete = complete, show.cases = show.cases, incl.cut1 = incl.cut1, 
                  incl.cut0 = incl.cut0, inf.test = inf.test)
    }
    
    initial.data <- data
    outcome.copy <- outcome
    
    if (grepl("[{]", outcome)) {
     
        outcome <- substr(outcome, 1, regexpr("[{]", outcome)[1] - 1)
        outcome.value <- as.numeric(substr(outcome.copy, regexpr("[{]", outcome.copy)[1] + 1, 
                                                         regexpr("}", outcome.copy)[1] - 1))
     
        data[, outcome] <- ifelse(data[, outcome] == outcome.value, 1, 0)
    }
    
    if (all(gsub("\\s", "", exo.facs, perl = TRUE) == "")) {
     
        exo.facs <- colnames(data)[-which(colnames(data) == outcome)]
    } 
    
    # make sure all relevant labels are uppercase
    colnames(data[, c(exo.facs, outcome)]) <- toupper(colnames(data[, c(exo.facs, outcome)]))
    exo.facs <- toupper(exo.facs)
    outcome <- toupper(outcome)
    
    data <- data[, c(exo.facs, outcome)]
    colnames(initial.data) <- toupper(colnames(initial.data))
    initial.data <- initial.data[, c(exo.facs, outcome)]
    
    # if user lowers only incl.cut1, make sure incl.cut0 is lowered accordingly
    if (incl.cut0 > incl.cut1) {
     incl.cut0 <- incl.cut1
    }
    
    # If 'neg.out = TRUE', us the negation of the outcome
    if (neg.out) {
        data[, outcome] <- 1 - data[, outcome]
    }
    
    dc.code <- unique(unlist(lapply(data, function(x) {
        if (is.numeric(x)) {
            return(x[x < 0])
        }
        else {
            return(as.character(x[x %in% c("-", "dc")]))
        }
    })))
    
    if (length(dc.code) == 0) {
        dc.code <- -1
    }
    else if (length(dc.code) > 1) {
        stop("\nMultiple \"don't care\" codes found.\n\n", call. = FALSE)
    }
    
    data <- as.data.frame(lapply(data, function(x) {
        x <- as.character(x)
        x[x == dc.code] <- -1
        return(as.numeric(x))
    }))
    
    names(data) <- c(exo.facs, outcome)
    
    data[data < 0] <- -1
    rownames(data) <- rownames(initial.data)
    
    nofexo.facs <- length(exo.facs)
    fuzzy.cc <- apply(data[, exo.facs, drop = FALSE], 2, function(x) any(x %% 1 > 0))
    
    for (i in seq(length(exo.facs))) {
        
         if (!fuzzy.cc[i]) {
             copy.cc <- data[, i]
           
               if (any(copy.cc < 0)) {
               
                   copy.cc[copy.cc < 0] <- max(copy.cc) + 1
                   data[, i] <- copy.cc
               }
         }
    }
    
    # levels must begin with 0 and MUST be incremented by 1
    # perhaps trying something like
    # apply(data[, exo.facs], 2, function(x) length(unique(x))) + 1
    noflevels <- apply(data[, exo.facs, drop = FALSE], 2, max) + 1
    noflevels[noflevels == 1] <- 2
    noflevels[fuzzy.cc] <- 2
    
    if (via.pof) {
        return(as.vector(noflevels))
    }
    
    if (memcare) {
        mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
        inclpri <- .Call("truthTableMem", as.matrix(data[, exo.facs]), noflevels, mbase, as.numeric(fuzzy.cc), data[, outcome], package="QCApro")
    }
    else {
        tt <- mintermMatrix(noflevels)
        inclpri <- .Call("truthTable", as.matrix(data[, exo.facs]), tt, as.numeric(fuzzy.cc), data[, outcome], package="QCApro")
    }
    
    colnames(inclpri[[1]]) <- seq_len(ncol(inclpri[[1]]))
    
    if ("SCTT" %in% names(other.args)) {
        copyinclpri <- inclpri
    }
    
    line.data <- inclpri[[2]]
    
    preserve <- inclpri[[1]][3, ] >= n.cut
    inclpri  <- inclpri[[1]][1:2, ]
    # inclpri  <- inclpri[[1]][1:2, ncut.fre >= n.cut]
    
    inclpri[is.na(inclpri)] <- NA
    
    outvalues <- as.numeric(inclpri[1, preserve] >= (incl.cut1 - .Machine$double.eps ^ 0.5))
    outvalues[inclpri[1, preserve] < incl.cut1 & inclpri[1, preserve] >= (incl.cut0 - .Machine$double.eps ^ 0.5)] <- "C"
    names(outvalues) <- colnames(inclpri)[preserve]
    
    freq.lines <- table(line.data)
    
    line.data[!line.data %in% colnames(inclpri)[preserve]] <- 0
    
    excluded <- line.data == 0
    line.data <- line.data[!excluded]
    
    if (memcare) {
        data[!excluded, exo.facs] <- getRow(noflevels, line.data)
    }
    else {
        data[!excluded, exo.facs] <- tt[line.data, ]
    }
    
    
    if (any(excluded)) {
        excluded.cases <- data[excluded, ]
        data <- data[!excluded, ]
    }
    
    data[, outcome] <- outvalues[match(line.data, names(outvalues))]
    
    if (complete) { # implicitly memcare is FALSE
        line.tt <- seq_len(dim(tt)[1])
    }
    else { # this doesn't implicitly imply memcare = TRUE, although it could be
        line.tt <- sort(unique(line.data))
        tt <- getRow(noflevels, line.tt)
    }
    
    tt <- as.data.frame(tt)
    rownames(tt) <- line.tt
    
    tt <- cbind(tt, "?")
    colnames(tt) <- c(exo.facs, "OUT")
    tt$OUT <- as.character(tt$OUT)
    
    tt[, "n"] <- 0
    lines.in.tt <- which(names(freq.lines) %in% rownames(tt))
    tt[names(freq.lines)[lines.in.tt], "n"] <- freq.lines[lines.in.tt]
    
    
    outcome.values <- sort(unique(data[, outcome]))
    
    for (i in seq(length(outcome.values))) {
        linesubset <- table(line.data[data[, outcome] == outcome.values[i]])
        tt[match(names(linesubset), line.tt), "OUT"] <- outcome.values[i]
    }
    
    
    tt <- cbind(tt, incl="-", PRI="-")
    tt$incl <- as.character(tt$incl)
    tt$PRI <- as.character(tt$PRI)
    inclpri <- inclpri[, which(colnames(inclpri) %in% rownames(tt)), drop = FALSE]
    
    tt[colnames(inclpri), "incl"] <- inclpri[1, ]
    tt[colnames(inclpri), "PRI"] <- inclpri[2, ]
    
    if (any(sort.by != "")) {
        sort.args <- c("incl", "n")
        if (!all(is.na(args.match <- match(sort.by, sort.args)))) {
            sort.args <- sort.args[args.match[!is.na(args.match)]]
            consorder <- order(tt[, sort.args[1]], decreasing=decreasing)
            if (length(sort.args) == 2) {
                consorder <- order(tt[, sort.args[1]], tt[, sort.args[2]], decreasing=decreasing)
            }
            tt <- tt[consorder, ]
            line.tt <- line.tt[consorder]
        }
    }
    
    #return(list(line.tt, initial.data, line.data))
    
    cases <- sapply(line.tt, function(x) {
        paste(rownames(data)[which(line.data == x)], collapse=",")
    })
    
    for (i in seq(length(exo.facs))) {
        if (!fuzzy.cc[i]) {
            if (any(initial.data[, i] == dc.code)) {
                tt[, i][tt[, i] == max(tt[, i])] <- dc.code
                data[, i][data[, i] == max(data[, i])] <- dc.code
                noflevels[i] <- noflevels[i] - 1
            }
        }
    }
    
    statistical.testing <- FALSE
    
    if (inf.test[1] == "binom") {
        
        statistical.testing <- TRUE
        
        if (length(inf.test) > 1) {
        
            alpha <- as.numeric(inf.test[2])
        }
        
        else {
         
            alpha <- 0.05
        }
        
        observed <- which(tt$OUT != "?")
        trials <- tt[observed, "n"]
        success <- trials * as.numeric(tt[observed, "incl"])
                
        tt <- cbind(tt, pval1 = "-", pval0 = "-")
        tt[, "pval1"] <- tt[, "pval0"] <- as.character(tt[, "pval1"])
        tt[observed, "OUT"] <- 0
        
        for (i in seq_along(observed)) {
             
             pval1 <- tt[observed[i], "pval1"] <- binom.test(success[i], trials[i], incl.cut1, "less")$p.value
             pval0 <- tt[observed[i], "pval0"] <- binom.test(success[i], trials[i], incl.cut0, "greater")$p.value
            
             if (pval1 > alpha) {
            
                 tt[i, "OUT"] <- 1
             }
            
             else if (pval1 < alpha & pval0 < alpha) {
            
                 tt[i, "OUT"] <- "C"
             }
        }
        tt[!observed, "pval1"] <- "-"
        tt[!observed, "pval0"] <- "-"
    }
    
    if (show.cases) {
        tt <- cbind(tt, cases)
    }
    
    x <- list(tt = tt, indexes = sort(unique(line.data)), 
              noflevels = as.vector(noflevels), initial.data = initial.data, 
              recoded.data = data, cases = cases, neg.out = neg.out,
              inf.test = statistical.testing, incl.cut1 = incl.cut1, 
              incl.cut0 = incl.cut0)
    
    if (any(excluded)) {
       x$excluded <- excluded.cases
    }
    
    x$tt$incl[is.na(x$tt$incl)] <- "-"
    x$tt$PRI[is.na(x$tt$PRI)] <- "-"
    
    if (use.letters & sum(nchar(colnames(data)[-ncol(data)])) != (ncol(data) - 1)) { # also verify if not already letters
        colnames(x$tt)[seq(nofexo.facs)] <- LETTERS[seq(nofexo.facs)]
    }
    
    x$outcome <- outcome.copy
    
    PRI <- FALSE
    if ("PRI" %in% names(other.args)) {
        if (is.logical(other.args$PRI)) {
            PRI <- other.args$PRI[1]
        }
    }
    
    x$PRI <- PRI
    x$call <- metacall
    
    x$opts$use.letters <- use.letters
    
    if ("SCTT" %in% names(other.args)) {
        x$SCTT <- copyinclpri
    }
    
    return(structure(x, class = "tt"))
}

