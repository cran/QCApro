\name{mintermMatrix}

\alias{mintermMatrix}

\title{Create Minterm Matrices}

\description{
This function creates minterm and implicant matrices. It is mainly used for 
internal and demonstration purposes.
}

\usage{
mintermMatrix(noflevels, logical = FALSE)
}

\arguments{
  \item{noflevels}{The number of levels for each exogenous factor.}
  \item{logical}{Logical, return the matrix in logical values (only bivalent data).}
}

\details{
Minterm matrices contain all unique and complete conjunctions that can 
be formed from all levels of \eqn{k} factors (Dusa and Thiem 2015). The total 
number of minterms \eqn{d}{d} is given by 
\eqn{d = \prod_{j = 1}^{k}{p_{j}}}{d = \prod_{j = 1}^{k}{p_{j}}}, where 
\eqn{p_{j}}{p_{j}} is the number of levels for exogenous factor \eqn{j} and 
\eqn{k} is the total number of exogenous factors. A minterm matrix is an essential 
part of a truth table.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Dusa, Adrian, and Alrik Thiem. 2015. \dQuote{Enhancing the Minimization of Boolean and Multivalue Output Functions with \emph{e}QMC.} \emph{Journal of Mathematical Sociology} \strong{39} (2):92-108. DOI: \href{http://dx.doi.org/10.1080/0022250X.2014.897949}{10.1080/0022250X.2014.897949}.
}

\seealso{\code{\link{implicantMatrix}}, \code{\link{truthTable}}}

\examples{
# a minterm matrix with three bivalent exogenous factors
noflevels <- rep(2, 3)
mintermMatrix(noflevels)

# with logical values
mintermMatrix(noflevels, logical = TRUE)
}

\keyword{functions}
