`findSubsets` <-
function(noflevels3k, row.no, maximum) {
    maximum <- ifelse(missing(maximum), prod(noflevels3k), maximum)
    # print(noflevels3k - 1)
    # print(rev(c(1, cumprod(rev(noflevels3k))))[-1])
    result <- lapply(row.no, function(x) {
        .Call("findSubsets", x, noflevels3k - 1, rev(c(1, cumprod(rev(noflevels3k))))[-1], maximum)
    })
    
    return(sort(unique(unlist(result))))
}

