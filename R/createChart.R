`createChart` <- 
function (input, copyinput, rows, cols) {
	input2 <- matrix(logical(length(input)), dim(input))
	input2[input > 0] <- TRUE
	
	result <- sapply(seq(nrow(input)), function(x) {
		apply(copyinput, 1, function(y) {
			all(input[x, input2[x,]] == y[input2[x,]])
        })
    })
    
    if (nrow(copyinput) == 1) {
        result <- matrix(result)
    }
    else {
        result <- t(result)
    }
    
    if (!missing(rows)) rownames(result) <- rows
    if (!missing(cols)) colnames(result) <- cols
    return(result)
}

