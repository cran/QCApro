`fuzzyor` <- function(...) {
    other.args <- list(...)
    if (is.matrix(other.args[[1]]) | is.data.frame(other.args[[1]])) {
        return(apply(other.args[[1]], 1, max))
    }
    else {
        return(pmax(...))
    }
}
