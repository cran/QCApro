#`prettyString` <-
#function(preamble, mystring, blanks="") {
#    blankchars <- nchar(preamble) + nchar(blanks)
#    paste(paste(strwrap(paste(preamble, mystring),
#                        prefix=paste(rep(" ", blankchars), collapse=""),
#                        width=floor(getOption("width")*0.95), initial=""),
#                collapse="\n"), "\n", sep="")
#}


prettyString <- function(string.vector, string.width = 80, repeat.space = 5, 
                         separator = ",", sufnec = "", outcome = "", 
                         cases = FALSE) {
    
  if (length(string.vector) == 1) {
  
      if (nchar(paste0(string.vector, " ", sufnec, " ", outcome)) >= string.width) {

          string.vector <- unlist(strsplit(string.vector, split = paste0(" \\", separator, " ")))
      }
  }
    
  string <- string.vector[1]
    
  if (length(string.vector) > 1) {
  
      startpoint <- 1
   
      for (j in seq(2, length(string.vector) + 1)) {
      
           if (j <= length(string.vector)) {
                
               if (nchar(paste(string.vector[seq(startpoint, j - ifelse(separator == ";", 1, 0))], collapse = paste(ifelse(separator == ";", "", " "), separator, " ", sep=""))) >= string.width) {
                    # if (nchar(paste(string.vector[seq(startpoint, j - ifelse(separator == ";", 1, 0))], collapse = paste(ifelse(separator == ";", "", " "), separator, sep=""))) >= string.width) {
                    #     string <- paste(paste(string, "\n", sep=""), 
                    #                     paste(rep(" ", repeat.space), collapse=""),
                    #                     separator, " ", string.vector[j], sep="")
                    # }
                    # else {
                   string <- paste0(paste(string, ifelse(separator == ";", "", " "), separator, "\n", sep = ""), 
                               paste(rep(" ", repeat.space), collapse=""),
                                 string.vector[j])
                    # }
                    
                   startpoint <- j
               }
               
               else {
                   
                   string <- paste0(string, ifelse(separator == ";", "", " "), separator, " ", string.vector[j])
               }
           }
           
           else {
               
               if (outcome != "") {
                    
                    # j is already bigger than the length of the string.vector
                   last.part <- paste(paste(string.vector[seq(startpoint, j - 1)], collapse = paste(ifelse(separator == ";", "", " "), separator, " ", sep="")), sep="")
                    
                   if (nchar(paste0(last.part, " ", sufnec, " ", outcome)) >= string.width) {
                    
                       string <- paste0(paste0(string, "\n"),
                                   paste(rep(" ", repeat.space), collapse = ""),
                                     sufnec, " ", outcome)
                   }
                   
                   else {
                       
                       string <- paste0(string, " ", sufnec, " ", outcome)
                   }
               }
           }
       }
  }
    
  else {
       
      if (outcome != "") {
        
          string <- paste0(string, " ", sufnec, " ", outcome)
      }
  }
  
  return(string)
}


