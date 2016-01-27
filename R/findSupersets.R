findSupersets <- function (noflevels3k, input.combs) {
    if (!is.matrix(input.combs)) {
        if (!is.vector(input.combs)) {
            cat("\n")
            stop("input.combs must be either an solution-space matrix or a vector of line numbers.\n\n",
                 call. = FALSE)
        }
        else {
            if (any(input.combs > prod(noflevels3k))) {
                cat("\n")
                stop(paste("Some line numbers do not belong in the solution-space for",
                           length(noflevels3k), "causal conditions.\n\n"), call. = FALSE)
            }
            input.combs <- getRow(noflevels3k, input.combs)
        }
    }
    mbase <- rev(c(1, cumprod(rev(noflevels3k))))[-1]
    allcombn <- t(mintermMatrix(rep(2, length(noflevels3k)))[-1, ])
    primes <- sort(unique(as.vector(apply(input.combs, 1, function(x) (x*mbase) %*% allcombn + 1))))
    if (primes[1] == 1) {
        return(primes[-1])
    }
    else {
        return(primes)
    }
}

