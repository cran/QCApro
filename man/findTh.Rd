\name{findTh}

\alias{findTh}

\title{Find Calibration Thresholds}

\description{
This function finds calibration thresholds for splitting base variables into the 
desired number of groups using cluster analysis.
}

\usage{
findTh(x, groups = 2, hclustm = "complete", distm = "euclidean")
}

\arguments{
  \item{x}{An interval or ratio-scaled base variable.}
  \item{groups}{A vector of integers with the desired number of groups.}
  \item{hclustm}{The agglomeration (clustering) method to be used.}
  \item{distm}{The distance measure to be used.}
}

\details{
For more details about argument \code{groups}, see \code{?cutree}. For more 
details about argument \code{hclustm}, see \code{?hclust}. For more details 
about argument \code{distm}, see \code{?dist}.
}

\value{
A numeric vector of suggested threshold(s) for dividing base variables 
into the desired number of groups.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\note{Default values from the \code{hclust} method and the \code{dist} method 
are used for both the distance measure \code{distm} and the clustering method 
\code{hclustm}.
}

\seealso{\code{\link{cutree}}, \code{\link{hclust}}, \code{\link{dist}}}

\examples{
# 15 random values between 1 and 100 
x <- sample(1:100, size = 15)

# split into two groups for csQCA
findTh(x)

# split into three groups for mvQCA
findTh(x, groups = 3)
}

\keyword{functions}

