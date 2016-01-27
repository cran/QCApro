implicIndep <- function (expression, n.samples = 1, size.sample = 100, corr = "0") {
 
  if (missing(expression)) {
   
      errmsg <- paste0("The expression is missing.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  } 
  
  expression <- gsub(" ", "", expression)
  
  antec <- substr(expression, 2, regexpr("<", expression)[1] - 1)
  
  out.f <- substr(expression, regexpr(">", expression)[1] + 1, regexpr(")", expression)[1] - 1)
  
  pos.irr.f <- as.vector(gregexpr("[+)]", expression)[[1]])
  
  pos.irr.f <- tail(pos.irr.f, 2)
  
  irr.f <- toupper(substr(expression, pos.irr.f[1] + 1, pos.irr.f[2] - 1))
  
  # check whether tautology has syntax (x + X); swap to (X + x) if yes
  pos.last.lower <- tail(as.vector(gregexpr("[a-z]", expression)[[1]]), 1)
  pos.last.upper <- tail(as.vector(gregexpr("[A-Z]", expression)[[1]]), 1)
  
  if (tail(pos.last.lower, 1) < tail(pos.last.upper, 1)) {
      substr(expression, pos.last.lower, pos.last.lower) <- irr.f
      substr(expression, pos.last.upper, pos.last.upper) <- tolower(irr.f)
  }
  
  # extract factors and rearrange conjunctions
  disj.list <- as.list(unlist(strsplit(antec, "[+]")))
  
  conj.list <- lapply(disj.list, function (x) {
  
      if (grepl("[*]", antec)) {
          sort(unlist(strsplit(x, "[*]")))
      }
   
      else if (!grepl("[0-9]", antec)) {
          strs <- substring(x, 1:nchar(x), 1:nchar(x))
          sort(strs[strs != ""])
      } 
    
      else {sst <- strsplit(x, "")[[1]]
          sort(paste0(sst[c(TRUE, FALSE)], sst[c(FALSE, TRUE)]))
      }
    }
  )
 
  if (any(sapply(conj.list, function (x) any(duplicated(toupper(x)))))) {
   
      errmsg <- paste0("Check the expression. At least one disjunct contains
                        contradiction or the same conjunct more than once.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  conj.sorted <- sort(unlist(lapply(conj.list, paste, collapse = "*")))
 
  # build unevaluated truth table
  f.names <- c(unique(toupper(unlist(conj.list))), out.f, irr.f)
  tt <- data.frame(mintermMatrix(rep(2, length(f.names))))
  dimnames(tt) <- list(as.character(seq(2^length(f.names))), f.names)
 
  # expand expression and create evaluated truth table
  antec.sorted <- paste(conj.sorted, sep = "", collapse = "+")
  
  conj.repl <- gsub(
    "[*]", ",", paste0("pmin(", unlist(strsplit(antec.sorted, "[+]")), ")")
  )
  
  subs <- paste0("pmax(", paste0(conj.repl, collapse = ","), ")", "==", out.f)
  
  taut <- sub("[+]", ",", paste0(", pmax", substr(
    expression, regexpr(irr.f, expression)[1] - 1, nchar(expression)))
  )
 
  coll <- paste0(gsub("([A-Z]+)", "tt$\\U\\1",
                     gsub("(\\b[a-oq-z][a-z0-9]*)", "1-\\U\\1",
                          paste0("pmin(", paste0(subs, taut), ")"), perl = TRUE),
                     perl = TRUE), "==TRUE"
  )
 
  tt <- tt[eval(parse(text = coll)), ]
  
  out.col <- match(out.f, names(tt))
  
  irr.col <- match(irr.f, names(tt))
  
  if (identical(mean(tt[, out.col]), 1)) {
   
      errmsg <- paste0("Check the expression. The endogenous factor is a constant.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  if (size.sample < nrow(tt)) {
   
      errmsg <- paste0("The sample size must be at least ", nrow(tt), ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  # begin simulation
  ded <- 2^(length(f.names) - 1)
  dat.list <- sol.list <- vector("list", n.samples)
  
  if (corr == "0") {
      for (i in seq(length(dat.list))) {
           set.seed(i)
           dat.list[[i]] <- rbind(tt, head(tt[sample.int(nrow(tt), size.sample,
                              replace = TRUE), ], -ded))[order(sample.int(size.sample, size.sample)), ]
           
           rownames(dat.list[[i]]) <- as.character(seq(size.sample))
           
           sol.list[[i]] <- eQMC(dat.list[[i]], outcome = out.f)$solution
      }
  } 
  else {
      
      if (corr == "+") {
          add.to <- tt[tt[, out.col] == tt[, irr.col], ]
      } 
      
      else {
          add.to <- tt[tt[, out.col] != tt[, irr.col], ]
      }
      
      for (i in seq_along(dat.list)) {
           set.seed(i)
           dat.list[[i]] <- rbind(tt, head(add.to[sample.int(nrow(add.to), size.sample,
             replace = TRUE), ], -ded))[order(sample.int(size.sample, size.sample)), ]
           
           rownames(dat.list[[i]]) <- as.character(seq(size.sample))
           
           sol.list[[i]] <- eQMC(dat.list[[i]], outcome = out.f)$solution
      }
  }
  
  # split up all models, sort conjuncts, then disjuncts, and merge again
  for (i in seq_along(sol.list)) {
       for (j in seq_along(sol.list[[i]])) {
            
            sol.list[[i]][[j]] <- lapply(sol.list[[i]][[j]], function (x) {
            
              paste0(sort(unlist(
                
                if (grepl("[*]", x)) {
                    strsplit(x, "[*]")
                }
                
                else {
                    strsplit(x, "")
                }
              )), collapse = "*")}
            )
       }
  }
  
  for (i in seq_along(sol.list)) {
       for (j in seq_along(sol.list[[i]])) {
            sol.list[[i]][[j]] <- gsub(
              "[*]", "", paste0(sort(unlist(sol.list[[i]][[j]], recursive = FALSE)), 
                                collapse = "+")
            )
       }
  }
  
  cor.list <- sapply(dat.list, function (x) round(cor(x[, out.col], x[, irr.col]), 4))
 
  # arrange returned object and test whether "irr.f" has been eliminated
  sols.obj <- list(tt = tt, dat.list = dat.list, sol.list = sol.list, cor.list = cor.list,
                  
    test = if ((any(sapply(sol.list, function (x) {
               grepl(irr.f, x, ignore.case = TRUE)})) == FALSE) &
               (gsub("[*]", "", antec.sorted) %in% unique(unlist(sol.list)) == FALSE)) {
               paste(irr.f, "has been eliminated (expression was not a causal structure).")
               } 
           else if (any(sapply(sol.list, function (x) {
                          grepl(irr.f, x, ignore.case = TRUE)}
                        )) == FALSE) {
               paste(irr.f, "has been eliminated.")
           } 
           
           else {
               
               paste(irr.f, "has not been eliminated.")
           }
  )
     return(sols.obj)
}