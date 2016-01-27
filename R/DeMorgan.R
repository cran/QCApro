DeMorgan <- function(expression, and.split = "", use.tilde = FALSE) {
    
    # print("bla")
    
    # TO DO: capture and error the usage of both "cD" and "D*E" in the same expression 
    
    
    # STRUCTURE of the big.list
    
    # level 1: split by separate components
        # "A + B(C + D*~E)" has two components "A" and "B(C + D*~E)"
    
    # level 2: split by brackets
        # "B(C + D*~E)" has "B" and "C + D*~E"
    
    # level 3: split by "+"
        # "C + D*~E" has "C" and "D*~E"
    
    # level 4: split by "*"
        # "D*~E" has "D" and "~E"
    
    # level 5: split by "~" (the result is only a vector, not a list)
        # "~E" has "~" and "E"
    
        
        # big.list <- splitMainComponents(expression)
        # big.list <- splitBrackets(big.list)
        # big.list <- removeSingleStars(big.list)
        # big.list <- splitPluses(big.list)
        # big.list <- splitStars(big.list)
        # big.list <- splitTildas(big.list)
        # big.list <- solveBrackets(big.list)
        # big.list <- simplifyList(big.list)
        
    if (class(expression) == "DeMorgan") {
        expression <- paste(expression[[1]][[2]], collapse = " + ")
    }
        
    splitMainComponents <- function(expression) {
        ind.char <- unlist(strsplit(expression, split=""))
        
        # remove all spaces (or white space)
        ind.char <- ind.char[ind.char != " "]
        
        if (grepl("\\(", expression)) {
            # split the string in individual characters
        
            open.brackets <- which(ind.char == "(")
            closed.brackets <- which(ind.char == ")")
            
            invalid <- ifelse(grepl("\\)", expression), length(open.brackets) != length(closed.brackets), FALSE)
            
            if (invalid) {
                cat("\n")
                stop("Invalid expression, open bracket \"(\" not closed with \")\".\n\n", call. = FALSE)
            }
            
            
            all.brackets <- sort(c(open.brackets, closed.brackets))
            
            if (length(all.brackets) > 2) {
                for (i in seq(3, length(all.brackets))) {
                    if (all.brackets[i] - all.brackets[i - 1] == 1) {
                        open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                        closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                    }
                    
                    if (all.brackets[i] - all.brackets[i - 1] == 2) {
                        if (ind.char[all.brackets[i] - 1] != "+") {
                            open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                            closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                        }
                    }
                }
            }
            
            for (i in seq(length(open.brackets))) {
                plus.signs <- which(ind.char == "+")
                last.plus.sign <- plus.signs[plus.signs < open.brackets[i]]
                if (length(last.plus.sign) > 0) {
                    open.brackets[i] <- max(last.plus.sign) + 1
                }
                else {
                    if (1 == 1) {
                        open.brackets[i] <- 1
                    }
                }
                next.plus.sign <- plus.signs[plus.signs > closed.brackets[i]]
                if(length(next.plus.sign) > 0) {
                    closed.brackets[i] <- min(next.plus.sign) - 1
                }
                else {
                    closed.brackets[i] <- length(ind.char)
                }
            }
                        
            # create an empty list with at least 3 times as many components as number of open brackets (just to make sure I have enough)
            big.list <- vector(mode="list", length = length(open.brackets) + 2)
            
            if (length(open.brackets) == 1) {
                # there is only one open bracket
                if (open.brackets > 1) {
                    # there's something before that open bracket
                    big.list[[1]] <- paste(ind.char[seq(1, open.brackets - 2)], collapse = "")
                }
                nep <- min(which(unlist(lapply(big.list, is.null))))
                big.list[[nep]] <- paste(ind.char[seq(open.brackets, closed.brackets)], collapse = "")
                if (closed.brackets < length(ind.char)) {
                    # there is something beyond the closed bracket
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(closed.brackets + 2, length(ind.char))], collapse = "")
                }
            }
            else {
                for (i in seq(length(open.brackets))) {
                    if (i == 1) {
                        # check if there's anything meaningful before the FIRST bracket
                        # i.e. containing a "+" sign, like "A + B(C + D)"
                        # before the first bracket is "A + B", but only B should be multiplied with "C + D"
                        
                        if (open.brackets[1] > 1) {
                            # there is something before the first bracket
                            big.list[[1]] <- paste(ind.char[seq(1, open.brackets[1] - 2)], collapse = "")
                        }
                        
                        nep <- min(which(unlist(lapply(big.list, is.null))))
                        big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                        
                    }
                    else {
                        nep <- min(which(unlist(lapply(big.list, is.null))))
                        big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                        
                        if (i == length(closed.brackets)) {
                            if (closed.brackets[i] < length(ind.char)) {
                                # there is something beyond the last closed bracket
                                nep <- min(which(unlist(lapply(big.list, is.null))))
                        
                                big.list[[nep]] <- paste(ind.char[seq(closed.brackets[i] + 2, length(ind.char))], collapse = "")
                                
                            }
                        }
                        
                    }
                }
            }
            
            nulls <- unlist(lapply(big.list, is.null))
            
            if (any(nulls)) {
                big.list <- big.list[-which(nulls)]
            }
            
        }
        else {
            big.list <- vector("list", length = 1)
            big.list[[1]] <- paste(ind.char, collapse="")
        }
        
        return(big.list)
    }
    
    
    #####
    # split each main component by separating brackets components
    splitBrackets <- function(big.list) {
        return(lapply(big.list, function(x) {
            as.list(unlist(strsplit(unlist(strsplit(x, split="\\(")), split="\\)")))
        }))
    }
    
    
    
    #####
    # remove individual components with single "*" signs 
    removeSingleStars <- function(big.list) {
        return(lapply(big.list, function(x) {
            single.components <- unlist(lapply(x, function(y) {
                return(y == "*")
            }))
            return(x[!single.components])
        }))
    }
    
    
    
    
    #####
    # split by "+"
    splitPluses <- function(big.list) {
        return(lapply(big.list, function(x) {
            lapply(x, function(y) {
                plus.split <- unlist(strsplit(y, "\\+"))
                return(as.list(plus.split[plus.split != ""]))
            })
        }))
    }
    
    
    
    #####
    # split by "*"
    splitStars <- function(big.list) {
        return(lapply(big.list, function(x) {
            lapply(x, function(y) {
                lapply(y, function(z) {
                    star.split <- unlist(strsplit(z, ifelse(and.split == "", "", paste0("\\", and.split))))
                    star.split <- star.split[star.split != ""]
                    if (and.split == "") {
                        tilda <- star.split == "~"
                        if (any(tilda)) {
                            tilda.pos <- which(tilda)
                            if (max(tilda.pos) == length(star.split)) {
                                cat("\n")
                                stop(paste0("Unusual expression \"", z, "\": terminated with a \"~\" sign?\n\n"), call. = FALSE)
                            }
                            star.split[tilda.pos + 1] <- paste0("~", star.split[tilda.pos + 1])
                            star.split <- star.split[-tilda.pos]
                        }
                    }
                    
                    return(as.list(star.split[star.split != ""]))
                })
            })
        }))
    }
    
    
    
    
    #####
    # split by "~"
    splitTildas <- function (big.list) {
        return(lapply(big.list, function(x) {
            lapply(x, function(y) {
                lapply(y, function(z) {
                    lapply(z, function(w) {
                        if (grepl("~", w)) {
                            wsplit <- unlist(strsplit(w, split=""))
                            if (max(which(wsplit == "~")) > 1) {
                                cat("\n")
                                stop(paste("Unusual expression: ", w, ". Perhaps you meant \"*~\"?\n\n", sep=""), call. = FALSE)
                            }
                            else {
                                return(c("~", sub("~", "", w)))
                            }
                        }
                        else {
                            return(w)
                        }
                    })
                })
            })
        }))
    }
    
    
    
    
    ######
    # determine if and which main components have brackets, and SOLVE them
    solveBrackets <- function(big.list) {
        bracket.comps <- which(unlist(lapply(big.list, length)) > 1)
        if (length(bracket.comps) > 0) {
            for (i in bracket.comps) {
                lengths <- unlist(lapply(big.list[[i]], length))
                indexes <- mintermMatrix(lengths) + 1
                ncol.ind <- ncol(indexes)
                i.list <- vector("list", length = nrow(indexes))
                
                for (j in seq(length(i.list))) {
                    i.list[[j]] <- vector("list", length = prod(dim(indexes)))
                    start.position <- 1
                    
                    for (k in seq(ncol.ind)) {
                        for (l in seq(length(big.list[[i]][[k]][[indexes[j, k]]]))) {
                            i.list[[j]][[start.position]] <- big.list[[i]][[k]][[indexes[j, k]]][[l]]
                            start.position <- start.position + 1
                        }
                    }
                    
                    if (start.position <= length(i.list[[j]])) {
                        i.list[[j]] <- i.list[[j]][- seq(start.position, length(i.list[[j]]))]
                    }
                }
                
                
                big.list[[i]] <- list(i.list)
            }
        }
        
        return(big.list)
    }
    
    
    
    
    simplifyList <- function(big.list) {
        lengths <- unlist(lapply(big.list, function(x) length(x[[1]])))
    
        big.list.copy <- vector("list", length = sum(lengths))
        
        start.position <- 1
        
        for (i in seq(length(big.list))) {
            for (j in seq(lengths[i])) {
                big.list.copy[[start.position]] <- big.list[[i]][[1]][[j]]
                start.position <- start.position + 1
            }
        }
        return(big.list.copy)
    }
    
    
    
    
    negateValues <- function(big.list, tilda = TRUE) {
        lapply(big.list, function(x) {
            lapply(x, function(y) {
                if (tilda) {
                    if (length(y) > 1) {
                        y <- toupper(y[2])
                    }
                    else {
                        if (use.tilde) {
                            y <- c("~", toupper(y))
                        }
                        else {
                            y <- tolower(y)
                        }
                    }
                }
                else {
                    if (y == toupper(y)) {
                        if (use.tilde) {
                            y <- c("~", toupper(y))
                        }
                        else {
                            y <- tolower(y)
                        }
                    }
                    else {
                        y <- toupper(y)
                    }
                }
            })
        })
    }
    
    
    
    
    removeDuplicates <- function(big.list) {
        
        big.list <- lapply(big.list, function(x) {
            
            values <- unlist(lapply(x, paste, collapse=""))
            x <- x[!duplicated(values)]

            
            # now trying to eliminate those which have both positive and negative
            # like "~A" and "A", or "a" and "A"
            ind.values <- unlist(x)
            ind.values <- ind.values[ind.values != "~"]
            ind.values <- toupper(ind.values)
            
            if (length(x) == 0 | any(table(ind.values) > 1)) {
                return(NULL)
            }
            else {
                return(x)
            }
        })
        
        big.list <- big.list[!unlist((lapply(big.list, is.null)))]
        
        
        # big.list.pasted
        blp <- lapply(big.list, function(x) {
            unlist(lapply(x, paste, collapse=""))
        })
        
        redundants <- vector(length = length(big.list))
        
        pairings <- combn(length(big.list), 2)
        
        for (i in seq(ncol(pairings))) {
            blp1 <- blp[[pairings[1, i]]]
            blp2 <- blp[[pairings[2, i]]]
            if (length(blp1) == length(blp2)) {
                if (all(sort(blp1) == sort(blp2))) {
                    redundants[pairings[2, i]] <- TRUE
                }
            }
            else {
                if (length(blp1) < length(blp2)) {
                    if (length(setdiff(blp1, blp2)) == 0) {
                        redundants[pairings[2, i]] <- TRUE
                    }
                }
                else {
                    if (length(setdiff(blp2, blp1)) == 0) {
                        redundants[pairings[1, i]] <- TRUE
                    }
                }
            }
        }
        
        return(big.list[!redundants])
        
    }
    
    
        
    if (is.qca(expression)) {
        result <- DeMorganLoop(expression)
    }
    else if (is.character(expression) & length(expression) == 1) {
        
        if (grepl("\\{", expression)) {
            if (grepl("~", expression)) {
                cat("\n")
                stop("Impossible combination of both \"~\" and \"{}\" multi-value notation.\n\n", call. = FALSE)
            }
            use.tilde <- FALSE
        }
        
        if (and.split == "" & grepl("\\*", expression)) {
            # cat("\n")
            # stop("The \"*\" symbol was found: consider using the argument and.split = \"*\".\n\n", call. = FALSE)
         and.split <- "*"
        }
        
        if (and.split != "" & and.split != "*") {
            if (!grepl(and.split, expression)) {
                cat("\n")
                stop("The AND-operator \"", and.split, "\" was not found.\n\n", call. = FALSE)
            }
        }
        
        big.list <- simplifyList(solveBrackets(splitTildas(splitStars(splitPluses(removeSingleStars(splitBrackets(splitMainComponents(expression))))))))
        
        flat.vector <- unlist(big.list)
        unique.values <- unique(flat.vector)
        
        already.letters <- all(nchar(unique.values) == 1)
        
        tilda <- ifelse(any(flat.vector == "~"), TRUE, FALSE)
        
        if (tilda) {
            use.tilde <- TRUE
        }
        
        if (tilda & and.split == "" & any(toupper(flat.vector) != flat.vector)) {
            cat("\n")
            stop("Unusual usage of both \"~\" sign and lower letters.\n\n", call. = FALSE)
        }
        
        negated.string <- paste("(", paste(unlist(lapply(negateValues(big.list, tilda), function(x) {
            paste(unlist(lapply(x, paste, collapse = "")), collapse = " + ")
        })), collapse = ")("), ")", sep="")
        
        
        big.list <- simplifyList(solveBrackets(splitTildas(splitStars(splitPluses(removeSingleStars(splitBrackets(splitMainComponents(negated.string))))))))
        
        
        # big.list <- splitMainComponents(negated.string)
        # big.list <- splitBrackets(big.list)
        # big.list <- removeSingleStars(big.list)
        # big.list <- splitPluses(big.list)
        # big.list <- splitStars(big.list)
        # big.list <- splitTildas(big.list)
        # big.list <- solveBrackets(big.list)
        # big.list <- simplifyList(big.list)
        
        
        initial <- expression
        negated <- unlist(lapply(removeDuplicates(big.list), function(x) {
            copyx <- unlist(lapply(x, function(y) {
                y <- y[y != "~"]
            }))
            x <- x[order(copyx)]
            paste(unlist(lapply(x, paste, collapse="")), collapse = and.split)
        }))
        
        result <- list(S1 = list(initial, negated))
        
    }
    
    return(structure(result, class = "DeMorgan"))
}





DeMorganLoop <- function(qca.object) {
 and.split <- qca.object$opts$collapse
    
    if ("i.sol" %in% names(qca.object)) {
        result <- vector("list", length=length(qca.object$i.sol))
        for (i in seq(length(qca.object$i.sol))) {
            names(result) <- names(qca.object$i.sol)
            result[[i]] <- lapply(qca.object$i.sol[[i]]$solution, paste, collapse = " + ")
            for (j in length(result[[i]])) {
                result[[i]][j] <- DeMorgan(result[[i]][[j]], and.split)
            }
        }
    }
    else {
        result <- lapply(lapply(qca.object$solution, paste, collapse = " + "), function(x) {
            DeMorgan(x, and.split)[[1]]
        })
        names(result) <- paste("S", seq(length(result)), sep="")
    }
    return(result)
}




