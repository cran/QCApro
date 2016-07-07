\name{implicantMatrix}

\alias{implicantMatrix}

\title{Create Implicant Matrices}

\description{
This function creates implicant matrices. An implicant matrix consists of all 
truth table minterms and their subsets, including the empty set.
}

\usage{
implicantMatrix(noflevels, raw = FALSE, arrange = FALSE)
}

\arguments{
  \item{noflevels}{The number of levels for each exogenous factor.}
  \item{raw}{Logical, return implicant matrix with indicator for elimination.}
  \item{arrange}{Logical, arrange for easier visual inspection.}
}

\details{
An implicant matrix consists of all minterms and their subsets, including the 
empty set (Dusa 2007, 2010; Thiem and Dusa 2015). The number of implicants \eqn{q}{q} is given by 
\eqn{q = \prod_{j = 1}^{k}{(p_{j} + 1)}}{q = \prod_{j = 1}^{k}{(p_{j} + 1)}}, 
where \eqn{p_{j}}{p_{j}} is the number of levels for factor 
\eqn{j}{j} and \eqn{k}{k} is the total number of exogenous factors.

If \code{raw = TRUE}, the indicator for elimination (\code{-1}) is used. 
}

\value{
A matrix with \eqn{\prod_{j = 1}^{k}{(p_{j} + 1)}}{\prod_{j = 1}^{k}{(p_{j} + 1)}} 
rows and \eqn{k}{k} columns.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Dusa, Adrian. 2007. \emph{Enhancing Quine-McCluskey}. COMPASSS: Working Paper 
2007-49. URL: \url{http://www.compasss.org/wpseries/Dusa2007b.pdf}.

Dusa, Adrian. 2010. \dQuote{A Mathematical Approach to the Boolean Minimization 
Problem.} \emph{Quality & Quantity} 44 (1):99-113.
DOI: \href{http://dx.doi.org/10.1007/s11135-008-9183-x}{10.1007/s11135-008-9183-x}.

Dusa, Adrian, and Alrik Thiem. 2015. \dQuote{Enhancing the Minimization of Boolean and Multivalue Output Functions with \emph{e}QMC.} \emph{Journal of Mathematical Sociology} \strong{39} (2):92-108. DOI: \href{http://dx.doi.org/10.1080/0022250X.2014.897949}{10.1080/0022250X.2014.897949}.
}

\seealso{
\code{\link{mintermMatrix}}
}

\examples{
# three exogenous factors with two levels each;
# first row is empty set
implicantMatrix(noflevels = rep(2, 3))

# two exogenous factors with three levels each
implicantMatrix(noflevels = rep(3, 2))

# arranged differently
implicantMatrix(noflevels = rep(3, 2), arrange = TRUE)

# with internal indicator for eliminated values
implicantMatrix(noflevels = rep(3, 2), raw = TRUE)
}

\keyword{functions}
