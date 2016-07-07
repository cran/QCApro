\name{randomDGS}
\alias{randomDGS}

\title{
Build a Random Data-Generating Structure
}

\description{
This evaluation function can be used to randomly build data-generating structures. It has initially been programmed for \href{https://www.researchgate.net/publication/280555425_Often_Trusted_But_Never_Properly_Tested_Evaluating_Qualitative_Comparative_Analysis}{Baumgartner and Thiem (2015)} to test the correctness of QCA's three search strategies (conservative/complex, intermediate, parsimonious).
}

\usage{
randomDGS(n.DGS = 1, exo.facs = c(""), seed.1 = NULL, seed.2 = NULL)
}

\arguments{
  \item{n.DGS}{The number of random data-generating structures to be built.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{seed.1}{The seed for the random generation of output function values.}
  \item{seed.2}{The seed for the random selection of a DGS in cases of structural
                ambiguities.}
}

\details{
The arguments .
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
Baumgartner, Michael, and Alrik Thiem. 2015. \emph{Often Trusted but Never 
(Properly) Tested: Evaluating Qualitative Comparative Analysis}. Paper presented 
at the 12th Conference of the European Sociological Association, 25-28 August, 
Czech Technical University, Prague (Czech Republic). \href{https://www.researchgate.net/publication/280555425_Often_Trusted_But_Never_Properly_Tested_Evaluating_Qualitative_Comparative_Analysis}{Link}.
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

