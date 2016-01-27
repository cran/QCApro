`findTh` <-
function(x, groups = 2, hclustm = "complete", distm = "euclidean") {
    x <- sort(x)
    cutpoints <- cumsum(rle(cutree(hclust(dist(x, method = distm), method = hclustm), k = groups))[[1]])
    values <- rep(NA, groups - 1)
    for (i in seq(length(values))) {
        values[i] <- mean(x[seq(cutpoints[i], cutpoints[i] + 1)])
    }
    return(values)
}

