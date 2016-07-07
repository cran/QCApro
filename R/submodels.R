submodels <- function (expression, noflevels = c(), test = TRUE) {
  
  function.call <- match.call()
  qca.object <- is.qca(expression)
  
  # test whether "expression" is an object of class 'qca'  
  if (qca.object) {
      
      cn <- colnames(expression$tt$tt)
            
      if (expression$opts$use.letters == TRUE | 
          all(cn[1:(match("OUT", cn) - 1)] %in% LETTERS)) {
      
          return(submodels.loop(expression = expression))
      }
   
      else {
       
          errmsg <- paste0("If an object of class 'qca' is passed to the 'submodels'
                            function, please ensure that all exogenous factors are 
                            labelled with single letters or that the argument 
                            'use.letters' is set to 'TRUE' in the call to the 
                            'eQMC' function.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      } 
  }

  expression.initial <- expression
 
  # eliminate all white space
  expression <- gsub(" ", "", expression)
  
  # extract antecedent if expression is minimally necessary
  if (grepl("<", expression)) {
  
      antec <- substr(expression, 1, regexpr("<", expression)[1] - 1)
  }
  
  # if it is not minimally necessary
  else if (grepl("-|=", expression)) {
    
      antec <- substr(expression, 1, regexpr("-|=", expression)[1] - 1)
  }
  
  # and if it already consists of the antecedent only
  else {
   
      antec <- expression
  }
  
  star <- grepl("[*]", antec)
  
  if (star & grepl("[a-zA-Z]{2,}", antec)) {
   
      errmsg <- paste0("The conjunction operator is not used consistently in
                        the expression ", expression.initial,". Either the 
                        '*'-operator must be put between each single-letter 
                        conjunct or must not be used at all.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }

  # Is a right-pointing arrowhead contained in the expression? If not,...
  if (!grepl(">", expression) & grepl("=|-|<", expression)) {
   
      errmsg <- paste0("The form of the expression was not as expected. Please 
                        consult the documentation of the 'submodels' function.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  # ...if yes, use the original outcome if provided
  else if (grepl(">", expression)) {
   
      out.f <- substr(expression, regexpr(">", expression)[1] + 1, 
                      nchar(expression)
      )
  }
  
  # otherwise, create a pseudo outcome
  else {
   
      out.f <- "PSEUDO"
  }
  
  if (length(noflevels) == 0) {
   
      f.labels <- sort(unique(toupper(unlist(strsplit(gsub("[*]|[+]", "", antec), "")))))
      noflevels <- rep(2, length(f.labels) + 1)
      names(noflevels) <- c(f.labels, toupper(out.f))
  }
  
  else {
   
      noflevels <- c(noflevels, PSEUDO = 2)
  }
  
  if (any(duplicated(names(noflevels)))) {
      
      f.duplic <- names(noflevels)[which(duplicated(names(noflevels)))]
      errmsg <- paste0("The same factor (", f.duplic, ") must not appear in the 
                        antecedent and the outcome.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  split.up <- function (antec, star) {
   
    # create a list of all disjuncts 
    disjuncts <- unlist(strsplit(antec, "[+]"))
   
    # create a list of all conjuncts depending on the structure of the antecedent;
    # the star operator may occur between single-digit factor levels or multi-value
    # terms if provided by the user as a string, but it must occur between 
    # non-multi-value, multi-digit factor levels; it will never occur for objects 
    # of class 'qca' with option "use.letters = TRUE" or when all factors have already
    # single letters  
    conjuncts.base <- lapply(disjuncts, function (x) {
    
      # if a star operator is used,...
      if (star) {
     
          sort(unlist(strsplit(x, "[*]")))
      }
    
      # if there is no star operator, factor levels are labelled with one letter
      else {
     
          strs <- substring(x, 1:nchar(x), 1:nchar(x))
          sort(strs[strs != ""])
      }
    })
    
    if (any(sapply(conjuncts.base, function (x) any(duplicated(x))))) {
     
        errmsg <- paste0("Check the expression. At least one disjunct contains 
                          the same conjunct more than once.")
        cat("\n")
        stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
    }
    
    else if (length(conjuncts.base) > 1 & any(sapply(conjuncts.base, function (x) {
               any(duplicated(toupper(x)))}))
             ) {
     
        errmsg <- paste0("Check the expression. At least one disjunct contains a 
                          conjunct and its negation.")
        cat("\n")
        stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
    }
    
    else {
     
        return(conjuncts.base)
    }
  }
  
  # following code snippet by http://stackoverflow.com/users/211116/spacedman;
  # augmented with own additions
  doChar <- Vectorize(function(c) {
   
    sprintf("tt$%s==%s", toupper(c), ifelse(c %in% LETTERS, "1", "0"))                      }
  )
  
  doWord <- Vectorize(function(W) {
   
    cs <- strsplit(W, "")[[1]]
    paste0("(", paste(doChar(cs), collapse = " & "), ")")
  })
  
  processString <- function(antec, out.f, star) {
   
    parts <- if (!star) {
    
      strsplit(antec, "\\+")[[1]] 
    }
   
    else {
    
      strsplit(gsub("[*]", "", antec), "\\+")[[1]]
    }
   
    paste0("(", paste0(doWord(parts), collapse = " | "), ") == (", 
           ifelse(all(strsplit(out.f, "")[[1]] %in% letters), "!", ""), "tt$", 
           toupper(out.f), ")")
  }
  
  # test whether 'expression' is a causal model to begin with
  if (test) {
   
      tt <- data.frame(mintermMatrix(noflevels))
      colnames(tt) <- names(noflevels)
    
      antec.sorted.test <- paste0(sort(unlist(lapply(
        (split.up(antec = antec, star = star)), paste0, collapse = ifelse(star, "*", "")
        ))), collapse = "+"
      )
  
      expr <- processString(antec = antec.sorted.test, out.f = out.f, star = star)
      
      tt <- tt[eval(parse(text = expr)), ]
  
      if ((all(strsplit(out.f, "")[[1]] %in% LETTERS) & all(tt[ , toupper(out.f)] == 1)) |
          (all(strsplit(out.f, "")[[1]] %in% letters) & all(tt[ , toupper(out.f)] == 0))) {
   
          errmsg <- paste0("The provided model, ", expression.initial, ", 
                            is no causal structure. The antecedent is a tautology.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
      
      else if ((all(strsplit(out.f, "")[[1]] %in% LETTERS) & all(tt[ , toupper(out.f)] == 0)) |
               (all(strsplit(out.f, "")[[1]] %in% letters) & all(tt[ , toupper(out.f)] == 1))) {
       
          errmsg <- paste0("The provided model, ", expression.initial, ", 
                            is no causal structure. The antecedent is a contradiction.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
      
      else {
       
          sol.list.test <- list()
       
          sol.list.test <- eQMC(tt, outcome = toupper(out.f), 
                                neg.out = ifelse(all(strsplit(out.f, "")[[1]] %in% letters), 
                                TRUE, FALSE))$solution
       
          for (j in seq_along(sol.list.test)) {
        
               sol.list.test[[j]] <- lapply(sol.list.test[[j]], function (z) {
         
                 paste0(sort(unlist(
          
                   if (grepl("[*]", z)) {
           
                       strsplit(z, "[*]")
                   }
          
                   else {
           
                       strsplit(z, "")
                   }
                 )), collapse = ifelse(star, "*", ""))
               })
          }
       
          for (j in seq_along(sol.list.test)) {
        
               sol.list.test[[j]] <- paste0(
               sort(unlist(sol.list.test[[j]], recursive = FALSE)), 
                 collapse = "+")
          }
       
          if (!antec.sorted.test %in% unlist(sol.list.test)) {
        
              errmsg <- paste0("The provided model, ", expression.initial, ", 
                                is no causal structure. It is quite possibly a 
                                case of F + fG = F + G, FG + fH + GH = FG + fH 
                                or FG + fg + HF + HG = g + H.")
              cat("\n")
              stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
          }
  
      }
  }    

  conjuncts.base <- split.up(antec = antec, star = star)
  
  # form all combinations of conjuncts within each list element
  nofchars <- sapply(conjuncts.base, length)
  combine <- function(a, b) combn(b, a, paste, collapse = "*")
  conjuncts.all <- conjuncts.base
  
  for (i in seq_along(nofchars)) {
    
       conjuncts.all[[i]] <- unlist(lapply(0:nofchars[i], combine, conjuncts.base[[i]]))
  }
  
  if (length(conjuncts.all) > 1) {
  
      # create the full list of submodels  
      mods <- t(apply(expand.grid(conjuncts.all, KEEP.OUT.ATTRS = FALSE,
                                  stringsAsFactors = FALSE), 1, c))
  
      # remove sub-models that contain duplicated disjuncts within model
      mods <- mods[!apply(mods, 1, function (x) {
        any(duplicated(x, incomparables = ""))}), ]
  
      # within rows, order along number of conjuncts, then within matrix along 
      # all columns
      mods <- t(apply(mods, 1, function (x) {x[order(nchar(x), x)]}))
      mods <- mods[do.call(order, as.list(as.data.frame(mods))), ]
      
      # remove duplicated submodels
      mods <- mods[!duplicated(mods), ]
      
      # remove submodels that contain disjuncts that are in a subset relation  
      cbs <- combn(ncol(mods), 2)
      
      # auxiliary truth matrix recording truth values for whether a subset 
      # relation between two disjuncts exists
      subs <- matrix(numeric(ncol(cbs)*nrow(mods)), nrow = nrow(mods))
  
      for (i in seq(nrow(subs))) {
    
        for (j in seq(ncol(cbs))) {
      
          if (mods[i, cbs[1, j]] != "") {
              subs[i, j] <- all(strsplit(mods[i, cbs[1, j]], split = "")[[1]] %in%
                                strsplit(mods[i, cbs[2, j]], split = "")[[1]])
          } 
      
          else {
       
              subs[i, j] <- 0
          }
        }
      }
      
      # only keep rows for which subs only has FALSE entries
      mods <- mods[rowSums(subs) < 1, ]
      
      # new auxiliary truth matrix recording truth values for submodels containing 
      # FG + fGH, including FG + fG, F + fG and F + f  
      subs <- matrix(numeric(ncol(cbs)*nrow(mods)), nrow = nrow(mods))
  
      for (i in seq(nrow(subs))) {
  
        for (j in seq(ncol(cbs))) {
  
          if (mods[i, cbs[1, j]] != "") {
          
              dis.a <- strsplit(mods[i, cbs[1, j]], split = "")[[1]]
              dis.b <- strsplit(mods[i, cbs[2, j]], split = "")[[1]]
              dis.ints <- intersect(dis.a, dis.b)
          
              if (length(setdiff(dis.a, dis.ints)) == 1) {
           
                  dis.a <- setdiff(dis.a, dis.ints)
                  dis.b <- setdiff(dis.b, dis.ints)
                  subs[i, j] <- tolower(dis.a) %in% tolower(dis.b)
              } 
          
              else {
             
                  subs[i, j] <- 0
              }
          } 
      
          else {
       
              subs[i, j] <- 0
          }
        }
      }

      mods <- mods[rowSums(subs) < 1, ]
      
      # as elimination by using syntactic rules with respect to disjuncts 
      # becomes too cumbersome, use 'eQMC' function to test whether a submodel 
      # is free of redundancies; first filter models with more than two disjuncts,
      # not all of which are also elementary conjuncts in the original antecedent
      # and at least one of which has a duplicated conjunct/negated conjunct
      idx <- which(rowSums(mods != "") > 2 & apply(mods, 1, function (x) {
        !all(x[x != ""] %in% unique(unlist(conjuncts.base))) &
           any(duplicated(toupper(strsplit(
             gsub("[*]","", paste(x, collapse = "")), "")[[1]]))
           )
        }
      ))

      mods.QMC.test <- t(apply(mods[idx, , drop = FALSE], 1, sort))
      
      if (length(idx) != 0) {
      
          mods <- mods[-idx, ]
      }
      
      subs <- numeric(nrow(mods.QMC.test))
  
      if (length(idx) != 0) {
      
          subs <- apply(mods.QMC.test, 1, function (x) {
            
            f.names <- c(unique(toupper(gsub("[*]", "", unlist(strsplit(x, "[*]"))))), 
                         toupper(out.f))
            
            tt <- data.frame(mintermMatrix(noflevels[match(f.names, names(noflevels))]))
            
            colnames(tt) <- f.names

            antec.eval <- paste(x[x != ""], collapse = "+")

            coll <- processString(antec = antec.eval, out.f = out.f, star = TRUE)

            tt <- tt[eval(parse(text = coll)), ]
     
            if (all(tt[, toupper(out.f)] == 1)) {
      
                return(2)
            }
     
            else {
       
                sol.list <- list()
            
                sol.list <- eQMC(tt, outcome = toupper(out.f), 
                                 neg.out = ifelse(all(strsplit(out.f, "")[[1]] %in% letters), 
                                   TRUE, FALSE))$solution
     
                for (j in seq_along(sol.list)) {
           
                     sol.list[[j]] <- lapply(sol.list[[j]], function (z) {
            
                       paste0(sort(unlist(
              
                         if (grepl("[*]", z)) {
                  
                             strsplit(z, "[*]")
                         }
              
                         else {
               
                             strsplit(z, "")
                         }
                       )), collapse = "*")
                     })
                }
     
                for (j in seq_along(sol.list)) {
                      
                     sol.list[[j]] <- paste0(
                       sort(unlist(sol.list[[j]], recursive = FALSE)), 
                       collapse = "+")
                }
     
                if (!antec.eval %in% unlist(sol.list)) {
           
                    return(1)
                } 
       
                else {
        
                    return(0)
                }
            }
          })
      }
  
      mods <- rbind(mods, mods.QMC.test[subs == 0, ])
      # redundancies <- mods.QMC.test[subs == 1, , drop = FALSE]
      # tautologies <- mods.QMC.test[subs == 2, , drop = FALSE]
  
      # flatten to vector  
      flatten <- function (x) {
        submods <- gsub(ifelse(star, "^[+]|[+]$", "^[+]|[*]|[+]$"), "",
                     gsub("[++]+", "+", apply(x, 1, function (y) {
                       paste(sort(y), collapse = "+")})))
        return(submods)
      }
  
      output.final <- list(expression = expression.initial, noflevels = noflevels, 
                           outcome = out.f, submodels = flatten(mods), 
                           # redundancies = flatten(redundancies),
                           # tautologies = flatten(tautologies), 
                           call = function.call)
  }
  
  else {
   
      output.final <- list(expression = expression.initial, noflevels = noflevels, 
                           outcome = out.f, 
                           submodels = gsub(ifelse(!star, "[*]", ""), "", 
                                            unlist(conjuncts.all)), 
                           # redundancies = character(0), tautologies = character(0),
                           call = function.call)
  }
  return(output.final)
}

 
submodels.loop <- function (expression) {
 
  model.list <- vector(mode = "list", length = length(expression$solution))
  lml <- length(model.list)
  names(model.list) <- paste0("M", seq(length(expression$solution)))

  incl.scores <- sapply(seq(lml), function (x) {
   
    expression$IC$individual[[x]]$sol.incl.cov[1]
  })
  
  cov.scores <- sapply(seq(lml), function (x) {
    
    expression$IC$individual[[x]]$sol.incl.cov[3]
  })
  
  relations <- sapply(seq(lml), function (x) {
   
    if (expression$relation == "suf" & cov.scores[x] < 1) {
   
        return(" => ")
    }
   
    else if (expression$relation == "suf" & cov.scores[x] == 1) {
    
        return(" <=> ")
    }
   
    else if (expression$relation == "sufnec" & cov.scores[x] >= incl.scores[x]) {
    
        return(" <=> ")
    }
    
    else {
    
      return(" => ")
    }
  })
  
  neg.out <- expression$tt$neg.out
  
  for (i in seq(length(expression$solution))) {

       antec.qca.model <- paste(expression$solution[[i]], collapse = " + ") 
       full.qca.model <- paste0(antec.qca.model, relations[i], 
                                ifelse(!neg.out, expression$tt$outcome, 
                                       tolower(expression$tt$outcome)))   
       
       # test = FALSE, otherwise, eQMC will throw an error message for antecedents 
       # consisting of single conditions testing the expression; in addition, no 
       # tests for objects of class 'qca' are needed because QCA models are directly 
       # causally interpretable
       model.list[[i]] <- submodels(expression = full.qca.model, test = FALSE)
  }
  
  return(model.list)
}
 
