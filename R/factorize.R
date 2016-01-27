factorize <- function(expression, and.split = "", sort.factorizing = FALSE, 
                      sort.factorized = FALSE) {
    
    factor.function <- function(trimmed.string) {
        
        my.string <- trimmed.string
        
        # create a list with all prime implicants split by literals
        if (and.split == "" & grepl("~", paste(trimmed.string, collapse = ""))) {
            
            list.my.string <- sapply(trimmed.string, strsplit, split = "")
            
            list.my.string <- lapply(list.my.string, function(x) {
                tildas <- x == "~"
                if (any(tildas)) {
                    x[which(tildas) + 1] <- paste("~", x[which(tildas) + 1], sep="")
                    x <- x[-which(tildas)]
                }
                return(x)
            })
        }
        else {
            list.my.string <- sapply(trimmed.string, strsplit, and.split)
        }
        
         # create a matrix with all combinations of prime implicants to be compared for similar literals
        all.combs <- mintermMatrix(rep(2, length(list.my.string)))
        all.combs <- all.combs[rowSums(all.combs) > 1, , drop=FALSE]
        all.combs <- col(all.combs) * as.vector(all.combs)
        
        
         # create a list with matched literals between prime implicants
        if (nrow(all.combs) > 1) {
         
            match.list <- as.list(apply(all.combs, 1, function(x) {
                x <- list.my.string[x[x > 0]]
                y <- table(unlist(x))
                return(names(y)[y == length(x)])
            }))
            names(match.list) <- lapply(match.list, paste, collapse = collapse)
        }
        else {
            match.list <- table(unlist(list.my.string))
            match.list <- list(names(match.list)[match.list == length(list.my.string)])
            names(match.list) <- lapply(match.list, paste, collapse = collapse)
        }
        
        if (length(match.list) > 0) {
             # see wich comparisons didn't yield similar literals
            null.branches <- unlist(lapply(match.list, function(x) all(is.na(x))))
             # erase those branches from the list
            match.list <- match.list[!null.branches]
            if (length(match.list) > 0) {
                
                if (nrow(all.combs) > 1) {
                     # and from all combinations
                    all.combs <- all.combs[!null.branches, , drop=FALSE]
                }
                
                if (sort.factorizing) {
                    sort.factorized <- FALSE
                    lengths.vector <- as.numeric(unlist(lapply(match.list, length)))
                    match.list <- match.list[rev(order(lengths.vector))]
                    all.combs <- all.combs[rev(order(lengths.vector)), ]
                }
                
                 # prepare a vector showing which columns from all.combs have been used
                 # to extract common.factor factors
                selected.rows <- rep(FALSE, nrow(all.combs))
                complex.list <- vector("list", length(selected.rows))
                
               
                extract <- function(match.list, all.combs, complex.list, my.string.index) {
                    initial.index <- my.string.index
                    for (i in 1:length(match.list)) {
                        common.factor <- match.list[[i]]
                         # see which other branches contain all common.factor literals from the current branch
                        similar.branches <- unlist(lapply(match.list[-i], function (x) all(common.factor %in% x)))
                        
                        if (any(similar.branches)) {
                             # see which are the other similar branches
                            similar.index <- seq(length(match.list))[-i][similar.branches]
                             # see which are the prime implicants with similar common.factor factors
                            my.string.index <- sort(unique(c(all.combs[c(i, similar.index), ])))
                            my.string.index <- my.string.index[my.string.index > 0]
                        }
                        else {
                             # see which are the prime implicants with similar common factors
                            my.string.index <- all.combs[i, ]
                            my.string.index <- my.string.index[my.string.index > 0]
                        }
                            
                         # paste the other literals from each index, separated by " + "
                        sol <- paste(sapply(my.string.index, function(x) {
                                paste(list.my.string[[x]][!list.my.string[[x]] %in% common.factor], collapse=collapse)
                                }), collapse=" + ")
                        
                        common.factor <- paste(match.list[[i]], collapse=collapse)
                        
                         # then combine everything having the common.factor in front of the paranthesys
                        factor.sol <- paste(common.factor, collapse, "(", sol, ")", sep="")
                        selected.rows <- apply(all.combs, 1, function(x) any(x %in% my.string.index))
                        
                        if (!is.null(initial.index)) my.string.index <- sort(unique(c(initial.index, my.string.index)))
                        
                        if (sum(!selected.rows) == 0) {
                             # no other comparison can be made; add all other prime implicants that have not been used
                            if (length(my.string[-my.string.index]) > 0) {
                                    factor.sol <- paste(factor.sol, paste(my.string[-my.string.index], collapse=" + "), sep=" + ")
                            }
                            names(complex.list)[i] <- factor.sol
                            complex.list[[i]] <- factor.sol
                        }
                        else {
                            sift <- function(x, y, z) {
                                sift.list <- list(match.list=NULL, all.combs=NULL)
                                sift.list[[1]] <- x[!z]
                                sift.list[[2]] <- y[which(!z), , drop=FALSE]
                                sift.list
                            }
                            sift.list <- sift(match.list, all.combs, selected.rows)
                            
                            names(complex.list)[i] <- factor.sol
                            complex.list[[i]] <- vector("list", length(sift.list$match.list))
                            complex.list[[i]] <- Recall(sift.list$match.list, sift.list$all.combs, complex.list[[i]], my.string.index)
                        }
                    }
                    return(complex.list)
                }
                
                
                my.string.index <- NULL
                complex.list <- extract(match.list, all.combs, complex.list, my.string.index)
                
                final.solution <- unique(names(unlist(complex.list)))
                
                if (length(final.solution) > 1) {
                    final.solution.list <- strsplit(final.solution, "\\.")
                    
                    if (sort.factorized) {
                        order.vector <- order(unlist(lapply(lapply(final.solution.list, "[", 1), nchar)), decreasing=TRUE)
                        final.solution.list <- final.solution.list[order.vector]
                        final.solution <- final.solution[order.vector]
                    }
                    
                    all.combs <- as.matrix(combn(length(final.solution.list), 2))
                    
                    match.list <- apply(all.combs, 2, function(x) {
                         # compare only solutions with the same length
                        if (length(final.solution.list[[x[1]]]) == length(final.solution.list[[x[2]]])) {
                             # return x (the indices from final.solution) if all items are equivalent
                            if (all(final.solution.list[[x[1]]] %in% final.solution.list[[x[2]]])) x
                        }
                    })
                    
                    # see if there are any null branches
                    null.branches <- unlist(lapply(match.list, is.null))
                    
                    if (!all(null.branches)) {
                         # remove those branches from match.list
                        match.list <- match.list[-which(null.branches)]
                         # the remaining branches contain equivalent (duplicated) solutions
                        equivalent.solutions <- unlist(lapply(match.list, "[", 2))
                         # remove equivalent solutions from final.solution
                        final.solution <- final.solution[-equivalent.solutions]
                    }
                    
                    final.solution <- gsub("\\.", " + ", final.solution)
                    
                }
                return(final.solution)
            }
        }
        else {
            return(NULL)
        }
    }
    
    
    getNonChars <- function(x) {
        # split by "+", incluging the trimming of the white space
        x <- gsub("^[[:space:]]+|[[:space:]]+$", "", unlist(strsplit(x, "\\+")))
        z <- vector(mode="list", length=length(x))
        for (i in seq(length(x))) {
            z[[i]] <- strsplit(gsub("[[:alnum:]]", "", x[i]), "+")[[1]]
        }
        z <- gsub("\\~", "", unique(unlist(z)))
        
        return(z[-which(z == "")])
    }
    
    collapse <- and.split
    if (and.split != "") and.split <- paste("\\", and.split, sep="")
    
    if (is.qca(expression)) {
        collapse <- and.split <- expression$opts$collapse
        if (and.split != "") and.split <- paste("\\", and.split, sep="")
        if ("i.sol" %in% names(expression)) {
            result <- list(i.sol=vector("list", length=length(expression$i.sol)))
            for (i in seq(length(expression$i.sol))) {
                names(result$i.sol) <- paste0(names(expression$i.sol), "S")
                result$i.sol[[i]] <- lapply(expression$i.sol[[i]]$solution, factor.function)
                names(result$i.sol[[i]]) <- unlist(lapply(expression$i.sol[[i]]$solution, paste, collapse=" + "))
            }
        }
        else {
            result <- lapply(expression$solution, function(x) {
                if (length(x) > 1) {
                    return(factor.function(x))
                }
                else {
                    return(NULL)
                }
            })
            names(result) <- unlist(lapply(expression$solution, paste, collapse=" + "))
        }
    }
    else if (is.DeMorgan(expression)) {
        
        if (names(expression)[1] == "S1") {
            result <- lapply(expression, function(x) {
                factor.function(x[[2]])
            })
            
            names(result) <- unlist(lapply(expression, function(x) {
                paste(x[[2]], collapse = " + ")
            }))
        }
        else {
            result <- list(lapply(expression, function(x) {
                int.result <- lapply(x, function(y) {
                    factor.function(y[[2]])
                })
                
                names(int.result) <- unlist(lapply(x, function(y) {
                    paste(y[[2]], collapse = " + ")
                }))
                
                return(int.result)
            }))
            names(result) <- "i.sol"
            names(result$i.sol) <- paste(names(result$i.sol), "N", sep="")
        }
        
    }
    else if (is.character(expression) & length(expression) == 1) {
        trimst <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)
        trimmed.str <- trimst(unlist(strsplit(expression, "\\+")))
        
        
        
        if (and.split != "") {
            if (!grepl(and.split, expression)) {
                cat("\n")
                stop("The product operator \"", and.split, "\" was not found.\n\n", call. = FALSE)
            }
        }
        else {
            nonchars <- getNonChars(trimmed.str)
            if (length(nonchars) > 0) {
                if (length(nonchars) > 1) {
                    cat("\n")
                    stop(paste("Multiple non alphanumeric characters found: \"", paste(nonchars, collapse=""), "\".\n\n", sep=""), call. = FALSE)
                }
                collapse <- nonchars
                and.split <- paste("\\", nonchars, sep="")
            }
        }
        
        if (length(trimmed.str) == 1) {
            result <- list(expression)
            names(result) <- expression
        }
        else {
            result <- list(factor.function(trimmed.str))
            names(result) <- expression
        }
    }
    
    return(structure(result, class = "fctr"))
}

