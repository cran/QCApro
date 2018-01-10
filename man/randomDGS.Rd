\name{randomDGS}
\alias{randomDGS}

\title{
Build a Random Data-Generating Structure
}

\description{
This evaluation function can be used to randomly build data-generating structures. It has initially been programmed for \href{http://dx.doi.org/10.1177/0049124117701487}{Baumgartner and Thiem (2017)} to test the correctness of QCA's three search strategies (conservative/complex, intermediate, parsimonious).
}

\usage{
randomDGS(n.DGS = 1, exo.facs = c(""), seed.1 = NULL, seed.2 = NULL, 
          prob = 0.5, diversity = 1, delete.trivial = FALSE)
}

\arguments{
  \item{n.DGS}{The number of random data-generating structures to be built.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{seed.1}{The seed for the random generation of output function values.}
  \item{seed.2}{The seed for the random selection of a DGS in cases of structural
                ambiguities.}
  \item{prob}{The probability of assigning a positive output function value to a minterm.}
  \item{diversity}{The diversity index value.}
  \item{delete.trivial}{Logical, delete "TRUE" and "FALSE" from set of structures.}
}

\details{
The argument \code{n.DGS} specifies the number of random data-generating structures to be built. 

The argument \code{exo.facs} is a character vector with the names of the exogenous factors. 

The argument \code{seed.1} sets the seed for the random generation of output function values, whereas \code{seed.2} sets the seed for the random selection of a DGS in cases of structural ambiguities. 

The argument \code{prob} is the probability of assigning a positive output function value to a minterm. 
The argument \code{diversity} specifies the diversity index value. It must be a number between 0 and 1.

The argument \code{delete.trivial} is logical, and specifies whether "TRUE" and "FALSE" should be deleted from the set of structures.
}

\value{
A list with the following two components:\cr
  \item{DGS}{A vector of the data-generating structure(s).}
  \item{tt}{The corresponding truth table(s).}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Thiem, Alrik        \tab: development, documentation, programming, testing
}}

\references{
Baumgartner, Michael, and Alrik Thiem. 2017. \dQuote{Often Trusted but Never 
(Properly) Tested: Evaluating Qualitative Comparative Analysis.} \emph{Sociological Methods & Research}. Advance online publication. DOI: \href{http://dx.doi.org/10.1177/0049124117701487}{10.1177/0049124117701487}.
}

\seealso{\code{\link{submodels}}}

\examples{
# randomly generate three data-generating structures on the basis of four
# exogenous factors
str <- randomDGS(n.DGS = 3, exo.facs = LETTERS[1:4], seed.1 = 1375, seed.2 = 3917)
str$DGS

# all correctness-preserving submodels of DGS 2, bd + abC, can then be found with the 
# 'submodels' function
submodels(str$DGS[2])$submodels
}

