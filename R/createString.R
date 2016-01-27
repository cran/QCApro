createString <- function(mydata, collapse = "*", uplow = FALSE, use.tilde = FALSE) {
    mydata <- changemydata <- as.matrix(mydata)
    exo.facs <- colnames(mydata)
    if (uplow) {
        changemydata[mydata == 0] <- tolower(rep(exo.facs, each=nrow(mydata))[mydata == 0])
        changemydata[mydata == 1] <- toupper(rep(exo.facs, each=nrow(mydata))[mydata == 1])
    }
    else if (use.tilde) {
        changemydata[mydata == 0] <- paste0("~", toupper(rep(exo.facs, each=nrow(mydata))[mydata == 0]))
        changemydata[mydata == 1] <- toupper(rep(exo.facs, each=nrow(mydata))[mydata == 1])
    }
    else {
        for (i in sort(unique(as.vector(mydata)))) {
            changemydata[mydata == i] <- paste0(rep(exo.facs, each=nrow(mydata))[mydata == i], "{", i, "}")
        }
    }
    
    input <- rep(NA, nrow(mydata))
    
    for (i in 1:nrow(mydata)) {
        input[i] <- paste(changemydata[i, ], collapse = collapse)
    }
    return(input)
}

