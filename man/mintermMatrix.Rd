\name{mintermMatrix}

\alias{mintermMatrix}

\title{Create Minterm and Implicant Matrices}

\description{
This function creates configuration and implicant matrices. It is mainly 
used for internal and demonstration purposes.
}

\usage{
mintermMatrix(noflevels, logical = FALSE)
}

\arguments{
  \item{noflevels}{The number of levels (values) for each condition variable.}
  \item{logical}{Logical, return the matrix in logical values (only bivalent data).}
}

\details{
Minterm matrices contain all unique and complete intersections that can 
be formed from all values of \eqn{k} variables. The total number of configurations
\eqn{d}{d} is given by 
\eqn{d = \prod_{j = 1}^{k}{p_{j}}}{d = \prod_{j = 1}^{k}{p_{j}}}, 
where \eqn{p_{j}}{p_{j}} is the number of values for condition variable 
\eqn{j} and \eqn{k} is the total number of condition variables. A 
configuration matrix is an essential part of a truth table.

An implicant matrix consists of all truth table configurations and their 
supersets, including the empty set (Dusa 2007, 2010). The number of implicants 
\eqn{q}{q} is given by 
\eqn{q = \prod_{j = 1}^{k}{(p_{j} + 1)}}{q = \prod_{j = 1}^{k}{(p_{j} + 1)}}.
For demonstration purposes, implicant matrices should be created by the 
function \code{implicantMatrix()}.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Dusa, Adrian. 2007. \emph{Enhancing Quine-McCluskey}. COMPASSS: Working Paper 
2007-49. URL: \url{http://www.compasss.org/wpseries/Dusa2007b.pdf}.

Dusa, Adrian. 2010. \dQuote{A Mathematical Approach to the Boolean Minimization 
Problem.} \emph{Quality & Quantity} 44 (1):99-113.
DOI: \href{http://dx.doi.org/10.1007/s11135-008-9183-x}{10.1007/s11135-008-9183-x}.
}

\seealso{\code{\link{implicantMatrix}}, \code{\link{truthTable}}}

\examples{
# a configuration matrix with three bivalent condition variables
noflevels <- rep(2, 3)
mintermMatrix(noflevels)

# with logical values
mintermMatrix(noflevels, logical = TRUE)

# its implicant matrix ("-1" stands for an eliminated value)
mintermMatrix(noflevels + 1) - 1

# a configuration matrix (the second variable has three values) 
noflevels <- c(2, 3, 2)
mintermMatrix(noflevels)

# its implicants matrix
mintermMatrix(noflevels + 1) - 1
}

\keyword{functions}
