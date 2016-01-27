# called by 'superSubset'
################################################################################

verify.data <- function(data, outcome = "", exo.facs = c("")) {
 
  # checking for absence of missing values
  if (any(is.na(data))) {
  
      nafactors <- names(data)[apply(apply(data, 2, is.na), 2, any)]
      errmsg <- paste0("The data must not contain missing values. The following
                   factors contain missing values: ", nafactors, ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
 } 
    
  # check the outcome specified by the user
  if (nchar(outcome) == 0) {
      
      stop("\nYou have not specified the outcome.\n\n")
  }
 
  else if (!outcome %in% colnames(data)) {
   
      errmsg <- paste0("The name of the outcome is incorrect.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
    
  # subset the data with the exo.facs specified
  if (length(exo.facs) > 1) {
      
      if (outcome %in% exo.facs) {
       
          errmsg <- paste0("The factor of the outcome ", outcome, " cannot also 
                            be an exogenous factor.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
   
      if (!all(exo.facs %in% names(data))) {
       
          errmsg <- paste0("The names of the exogenous factors are incorrect.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
  }
 
  else if (nchar(exo.facs) > 0) {
      
      if (outcome %in% exo.facs) {
       
          errmsg <- paste0("The factor of the outcome ", outcome, " cannot also 
                            be an exogenous factor.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
      
      else {
       
          errmsg <- paste0("At least two exogenous factors need to be specified.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
  }
}

# called by 'truthTable' and 'eQMC'
################################################################################

verify.tt <- function(data, outcome = "", exo.facs = c(""), complete = FALSE, 
                      show.cases = FALSE, incl.cut1 = 1, incl.cut0 = 1, 
                      inf.test = c(""), use.letters = FALSE) {
  
  # 'outcome' 
  #-----------------------------------------------------------------------------
  
  outcome.copy <- outcome
  # if the outcome is not specified,...
  if (gsub("\\s", "", outcome, perl = TRUE) == "") {
  
      errmsg <- paste0("No outcome is specified.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  # ...if it is, but the (bivalent) endogenous factor is not in 'data',...
  else if (!grepl("[{]", outcome) & !(outcome %in% colnames(data))) {
  
      errmsg <- paste0("The name of the outcome is incorrect. The factor ",
                        outcome, " does not exist in the data.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  # ...if it is, but the multivalent endogenous factor is not in 'data' or an
  # an incorrect level is specified,...
  else if (grepl("[{]", outcome)) {
  
      outcome <- substr(outcome, 1, regexpr("[{]", outcome)[1] - 1)
      outcome.value <- as.numeric(substr(outcome.copy, regexpr("[{]", outcome.copy)[1] + 1, 
                                                       regexpr("}", outcome.copy)[1] - 1))
  
      if (!(outcome %in% colnames(data))) {
   
          errmsg <- paste0("The endogenous factor ", outcome, " does not exist 
                            in the data.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), 
               call. = FALSE)
      }
  
      else if (!(outcome.value %in% unique(data[, outcome]))) {
   
          errmsg <- paste0("The endogenous factor ", outcome, " has no level {", 
                           outcome.value, "}.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), 
          call. = FALSE)
      }
  }
 
  # 'exo.facs' 
  #----------------------------------------------------------------------------- 
  # if no exogenous factors are specified, use all factors except that of the outcome
  if (all(gsub("\\s", "", exo.facs, perl = TRUE) == "")) {
  
      exo.facs <- colnames(data)[-which(colnames(data) == outcome)]
      
      if (length(exo.facs) > 16) {
       
          wrnmsg <- paste0("You are about to include ", length(exo.facs), " 
                            exogenous factors in the analysis. Minimization may
                            not be possible.")
          cat("\n")
          warning(paste(strwrap(wrnmsg, exdent = 7), collapse = "\n"), 
                  immediate. = TRUE, call. = FALSE)
      }
  } 
  
  # if there are at least two exogenous factors,...
  else if (length(exo.facs) > 1) {
  
      # and the endogenous factor is also in the set of exogenous factors,...
      if (outcome %in% exo.facs) {
   
          errmsg <- paste0("The factor of the outcome ", outcome, " cannot also 
                            be an exogenous factor.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
   
      # if at least one exogenous factor is not in the data,...
      if (!all(exo.facs %in% colnames(data))) {
    
          f.notindata <- exo.facs[which(!(exo.facs %in% colnames(data)))]
          errmsg <- paste0("The following exogenous factors are not present in 
                            the data: ", f.notindata, ".")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      } 
  }
 
  # if there is only one exogenous factor,...
  else if (nchar(exo.facs) > 0) {
  
      errmsg <- paste0("At least two exogenous factors need to be specified.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  data <- data[, c(exo.facs, outcome)]
  
  # if there are more than 26 exogenous factors (plus one outcome), letters
  # cannot be used
  if (use.letters & ncol(data) > 27) {
   
      errmsg <- paste0("Letters cannot be used. There are more than 26 
                       exogenous factors.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
 
  # missing data
  #-----------------------------------------------------------------------------
  # if there are missing values in the data under the factor frame,...
  if (any(is.na(data))) {
  
      f.nas <- names(data)[apply(apply(data, 2, is.na), 2, any)]
      errmsg <- paste0("Included factors must not contain missing values. The 
                        following factors contain missing values: ", f.nas, ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  } 
 
  # uncalibrated data
  #----------------------------------------------------------------------------- 
  # if the data under the factor frame has don't care values '-' or 'dc',
  # replace them with '-1'
  data <- as.data.frame(lapply(data, function(x) {
    
      x <- as.character(x)
      x[x %in% c("-", "dc")] <- -1
      return(as.numeric(x))
  }))
 
  # do the data contain values below -1 or values that have a modulo above 0
  # if they are larger than 0?
  f.uncalibrated <- names(data[ , sapply(data, function(x) {any(x < -1) | 
    (any(x %% 1 > 0 & x > 1))})]
  )
  
  # if there is at least one such factor found in the data,...
  if (!length(f.uncalibrated) == 0) {
  
      errmsg <- paste0("Uncalibrated data have been found in the following 
                        factors: ", paste(f.uncalibrated, collapse = " and "), ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
  
  # 'incl.cut1' and 'incl.cut0'
  #-----------------------------------------------------------------------------  
  # if one of the two inclusion cut-offs is below 0 or above 1,...
  if (any(c(incl.cut1, incl.cut0) < 0) | any(c(incl.cut1, incl.cut0) > 1)) {
      
      if ((incl.cut1 < 0 | incl.cut1 > 1) & (incl.cut0 >= 0 & incl.cut0 <= 1)) {
          
          errmsg <- paste0("The argument 'incl.cut1' has to be between 0 and 1. 
                            It is currently set to ", incl.cut1, ".") 
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
   
      else if ((incl.cut0 <  0 | incl.cut0 >  1) & 
               (incl.cut1 >= 0 & incl.cut1 <= 1)) {
          
          errmsg <- paste0("The argument 'incl.cut0' has to be between 0 and 1. 
                            It is currently set to ", incl.cut0, ".") 
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
   
      else {
   
          errmsg <- paste0("The arguments 'incl.cut1' and 'incl.cut0' have to be 
                            between 0 and 1. They are currently set to ", 
                           incl.cut1, " and ", incl.cut0, ", respectively.")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
  }
  
  # if incl.cut0 is EXPLICITLY set above incl.cut1,... 
  if (incl.cut0 < 1 & incl.cut0 > incl.cut1) {
      
      errmsg <- paste0("The value of the argument 'incl.cut0' must not be greater 
                        than that of the argument 'incl.cut1'. The former is currently 
                        set to ", incl.cut0, ", the latter to ", incl.cut1, ".")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }

  # if incl.cut1 is set below 0.5,... 
  if (incl.cut1 < 0.5) {
   
      wrnmsg <- paste0("The minimum sufficiency inclusion score for an output 
                        function value of '1' specified in the argument 
                        'incl.cut1' should not be below 0.5. It is currently set 
                        to ", incl.cut1, ".")
      cat("\n")
      warning(paste(strwrap(wrnmsg, exdent = 7), collapse = "\n"), 
              immediate. = TRUE, call. = FALSE)
  }
  
  # if incl.cut0 is set above 0.5,... 
  if (incl.cut0 > 0.5 & incl.cut1 > incl.cut0) {
   
      wrnmsg <- paste0("The maximum sufficiency inclusion score for an output 
                        function value of '0' specified in the argument 
                        'incl.cut0' should not be above 0.5. It is currently set 
                        to ", incl.cut0, ".")
      cat("\n")
      warning(paste(strwrap(wrnmsg, exdent = 7), collapse = "\n"), 
              immediate. = TRUE, call. = FALSE)
  }
 
  #-----------------------------------------------------------------------------
  # run tests for inf.test (see below)
  verify.inf.test(inf.test, data)
}

# called by 'eQMC'
################################################################################

verify.qca <- function(minimize = c("")) {
 
  # check if the user specifies something to minimize
  if (all(minimize == c(""))) {
    
      errmsg <- paste0("You have not specified any min-terms to be covered.
                        Normally, the setting is minimize = c(\"1\").")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
    
  # check if the user specifies something to minimize
  if (!all(minimize %in% c(0, 1, "C")) | all(c(0, 1) %in% minimize)) {

      errmsg <- paste0("The specified value to the argument 'minimize', ", 
                       minimize, ", is not allowed. Normally, the setting is 
                       minimize = c(\"1\").Please consult the documentation of 
                       the 'eQMC' function.")
      cat("\n")
      stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
  }
}

# called by 'eQMC'
################################################################################

verify.dir.exp <- function (data, outcome, exo.facs, dir.exp = c()) {
    
  # checking the directional expectations
  if (is.null(dir.exp)) {
  
       return(dir.exp)
  }
  
  else {
        
      # delc is dir.exp.list.complete
      delc <- vector("list", length = length(exo.facs))
      names(delc) <- exo.facs
        
      for (i in seq(length(delc))) {
            # sometimes a condition can have 0, 1 and "-" as values
            # see RagStr$EBA, which is also treated as a factor by default, in R
            
            values <- sort(unique(data[, exo.facs[i]]))
            if (is.factor(values)) {
                values <- as.character(values)
            }
            
            max.value <- values[length(values)]
            
            if (max.value != "-") {
                delc[[i]] <- rep(0, length(seq(0, as.numeric(max.value))))
                names(delc[[i]]) <- seq(0, as.numeric(max.value))
            }
        }
        
        if (length(dir.exp) != length(exo.facs)) {
            cat("\n")
            stop("Number of expectations does not match the number of exogenous factors.\n\n")
        }
        
        # del is dir.exp.list
        del <- strsplit(as.character(dir.exp), split=";")
        
        if (!is.null(names(dir.exp))) {
            if (length(names(dir.exp)) != length(exo.facs)) {
                cat("\n")
                stop("All directional expectations should have names, or none at all.\n\n")
            }
            else if (length(setdiff(names(dir.exp), exo.facs)) > 0) {
                cat("\n")
                stop("Incorrect names of the directional expectations.\n\n")
            }
            names(del) <- names(dir.exp)
            del <- del[exo.facs]
        }
        else {
            names(del) <- exo.facs
        }
        
        for (i in seq(length(del))) {
            
            values <- strsplit(del[[i]], split = "")
            values <- unlist(lapply(values, function(x) {
                paste(x[x != " "], collapse = "")
            }))
            
            if (all(values %in% c("-", "dc"))) {
                delc[[i]] <- -1
            }
            else {
                values <- setdiff(values, c("-", "dc"))
                if (length(setdiff(values, names(delc[[i]])) > 0)) {
                    cat("\n")
                    errmessage <- paste("Values specified in the directional expectations do not appear in the data, for condition \"", exo.facs[i], "\".\n\n", sep="")
                    stop(paste0(strwrap(errmessage, exdent = 7), collapse = "\n"))
                }
                else {
                    delc[[i]][as.character(values)] <- 1
                }
            }
        }
        return(delc)
    }
}

# called by 'eQMCLoop'
################################################################################

verify.mqca <- function(data, outcome = c(""), exo.facs = c("")) {
    
    mvoutcome <- grepl("[{]", outcome)
    outcome.value <- rep(-1, length(outcome))
    
    if (any(mvoutcome)) {
        outcome.copy <- outcome
        
        outcome.copy <- strsplit(outcome.copy, split = "")
        outcome.name <- outcome.value <- vector(mode="list", length = length(outcome))
        
        for (i in seq(length(outcome.copy))) {
            if (mvoutcome[i]) {
                outcome.value[[i]] <- as.numeric(outcome.copy[[i]][which(outcome.copy[[i]] == "{") + 1])
                outcome.name[[i]] <- paste(outcome.copy[[i]][seq(1, which(outcome.copy[[i]] == "{") - 1)], collapse="")
            }
            else {
                outcome.value[[i]] <- -1
                outcome.name[[i]] <- outcome[i]
            }
        }
        
        outcome <- unlist(outcome.name)
        
        if (length(intersect(outcome, names(data))) < length(outcome)) {
            outcome <- setdiff(outcome, names(data))
            cat("\n")
            errmessage <- paste("Outcome(s) not present in the data: \"", paste(outcome, collapse="\", \""), "\".\n\n", sep="")
            stop(paste0(strwrap(errmessage, exdent = 7), collapse = "\n"))
        }
        
        for (i in seq(length(outcome))) {
            if (mvoutcome[i]) {
                if (!any(unique(data[, outcome.name[[i]]]) == outcome.value[[i]])) {
                    cat("\n")
                    errmessage <- paste("The value {", outcome.value[[i]], "} does not exist in the outcome \"", outcome.name[[i]], "\".\n\n", sep="")
                    stop(paste0(strwrap(errmessage, exdent = 7), collapse = "\n"))
                }
            }
        }
        
        outcome.value <- unlist(outcome.value)
    }
    else {
        
        if (length(intersect(outcome, names(data))) < length(outcome)) {
            outcome <- setdiff(outcome, names(data))
            cat("\n")
            errmessage <- paste("Outcome(s) not present in the data: \"", paste(outcome, collapse="\", \""), "\".\n\n", sep="")
            stop(paste0(strwrap(errmessage, exdent = 7), collapse = "\n"))
        }
        
        fuzzy.outcome <- apply(data[, outcome, drop = FALSE], 2, function(x) any(x %% 1 > 0))
        
        # Test if outcomes are multivalent, even if the user did not specify this
        if (any(!fuzzy.outcome)) {
            outcome.copy <- outcome[!fuzzy.outcome]
            
            for (i in outcome.copy) {
                valents <- unique(data[, i])
                if (!all(valents %in% c(0, 1))) {
                    
                    errmessage <- paste0("Please specify the level of the endogenous 
                                          factor \"", i, "\" to used as the outcome .\n\n")
                    cat("\n")
                    stop(paste0(strwrap(errmessage, exdent = 7), collapse = "\n"))
                }
            }
        }
        
    }
    
    if (all(exo.facs == c(""))) {
        exo.facs <- colnames(data)
    }
    
    if (length(setdiff(outcome, exo.facs)) > 0) {
        outcome <- setdiff(outcome, exo.facs)
        cat("\n")
        errmessage <- paste("Outcome(s) not present in the conditions' names: \"", paste(outcome, collapse="\", \""), "\".\n\n", sep="")
        stop(paste(strwrap(errmessage, exdent = 7), collapse = "\n", sep=""))
    }
    
    invisible(return(list(mvoutcome = mvoutcome, outcome = outcome, outcome.value = outcome.value, exo.facs = exo.facs)))
    
}

# called by 'truthTable', 'eQMC'
################################################################################

verify.inf.test <- function(inf.test, data) {
 
  # is binomial testing specified?  
  if (all(inf.test != "")) {
      
      if (inf.test[1] != "binom") {
          
          errmsg <- paste0("Please check the first value to the argument 'inf.test'. 
                            It is currently set to ", inf.test[1], ". Only binomial 
                            testing is presently supported ('binom').")
          cat("\n")
          stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), call. = FALSE)
      }
    
      # if YES, do the data contain values 0 < x < 1?   
      else { 
     
          fuzzy <- apply(data, 2, function(x) any(x %% 1 > 0))
     
            if (any(fuzzy)) {
                
                errmsg <- paste0("The binomial test is unsuitable for fuzzy-set 
                                  data.")
                cat("\n")
                stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), 
                     call. = FALSE)
            }
      }
    
      # are both arguments to inf.test specified?  
      if (length(inf.test) > 1) {
       
          alpha <- as.numeric(inf.test[2])
         
          if (alpha < 0 | alpha > 1) {
                
              errmsg <- paste0("The critical significance level of the binomial 
                                test specified as the second value to the argument 
                                'inf.test' must be a number between 0 and 1. It 
                                is currently set to ", alpha, ".")
              cat("\n")
              stop(paste(strwrap(errmsg, exdent = 7), collapse = "\n"), 
                   call. = FALSE)
          }
      }
  }
}