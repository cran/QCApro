\name{testTESA}
\alias{testTESA}

\title{
Testing Schneider and Wagemann's Theory-Guided/Enhanced Standard Analysis (T/ESA)
}

\description{
This evaluation function can be used to test the implications of Schneider and Wagemann's 
Theory-Guided/Enhanced Standard Analysis (T/ESA; \href{http://dx.doi.org/10.1177/1065912912468269}{Schneider and Wagemann 2013}), and in particular, the procedure's first two stages, with respect to the extent of remainders that would have to be declared insufficient for the outcome. It has been programmed for \href{http://dx.doi.org/10.1093/pan/mpw024}{Thiem (2016)}.
}

\usage{
testTESA(data, outcome = "", neg.out = FALSE, exo.facs = c(""), n.cut = 1, 
         incl.cut1 = 1, incl.cut0 = 1)
}

\arguments{
  \item{data}{A dataset of bivalent crisp-set factors or bivalent fuzzy-set factors 
        or multivalent crisp-set factors.}
  \item{outcome}{The name of the outcome.}
  \item{neg.out}{Logical, use negation of \code{outcome}.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{n.cut}{The minimum number of cases with set membership score above 0.5 
        for an output function value of "0", "1" or "C".}
  \item{incl.cut1}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{incl.cut0}{The maximum sufficiency inclusion score for an output function 
        value of "0".}
}

\details{
The arguments \code{data}, \code{outcome}, \code{exo.facs}, \code{n.cut}, \code{incl.cut1} and \code{incl.cut0} are those of the \code{\link{eQMC}} function.
}

\value{
A numeric vector with the percentages of remainder minterms that would have been used as simplifying assumptions by Quine-McCluskey optimization but that were declared to be insufficient for the outcome by T/ESA.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Thiem, Alrik        \tab: development, documentation, programming, testing
}}

\references{
Ragin, Charles C. 2009. \dQuote{Qualitative Comparative Analysis Using Fuzzy Sets (fsQCA).} In \emph{Configurational Comparative Methods: Qualitative Comparative Analysis (QCA) and Related Techniques}, ed. B. Rihoux and C. C. Ragin. London: Sage Publications, pp. 87-121.

Schneider, Carsten Q., and Claudius Wagemann. 2013. \dQuote{Doing Justice to 
Logical Remainders in QCA: Moving Beyond the Standard Analysis.} \emph{Political 
Research Quarterly} \strong{66} (1):211-20. DOI: \href{http://dx.doi.org/10.1177/1065912912468269}{10.1177/1065912912468269}.

Thiem, Alrik. 2016. \dQuote{Standards of Good Practice and the Methodology of Necessary Conditions in Qualitative Comparative Analysis.} \emph{Political Analysis} \strong{24} (4):478-84. DOI: \href{http://dx.doi.org/10.1093/pan/mpw024}{10.1093/pan/mpw024}.  
}

\seealso{\code{\link{eQMC}}, \code{\link{truthTable}}, \code{\link{superSubset}}}

\examples{
# Schneider and Wagemann (2013, 212), using data from Ragin 
# (2009, 95), only present L and S as minimally necessary conditions
#-------------------------------------------------------------------

LIP <- data.frame(
 D = c(0.81,0.99,0.58,0.16,0.58,0.98,0.89,0.04,0.07,
       0.72,0.34,0.98,0.02,0.01,0.01,0.03,0.95,0.98),
 U = c(0.12,0.89,0.98,0.07,0.03,0.03,0.79,0.09,0.16,
       0.05,0.10,1.00,0.17,0.02,0.03,0.30,0.13,0.99),
 L = c(0.99,0.98,0.98,0.98,0.99,0.99,0.99,0.13,0.88,
       0.98,0.41,0.99,0.59,0.01,0.17,0.09,0.99,0.99),
 I = c(0.73,1.00,0.90,0.01,0.08,0.81,0.96,0.36,0.07,
       0.01,0.47,0.94,0.00,0.11,0.00,0.21,0.67,1.00),
 G = c(0.43,0.98,0.91,0.91,0.58,0.95,0.31,0.43,0.13,
       0.95,0.58,0.99,0.00,0.01,0.84,0.20,0.91,0.98),
 S = c(0.05,0.95,0.89,0.12,0.77,0.95,0.05,0.06,0.42,
       0.92,0.05,0.95,0.12,0.05,0.21,0.06,0.95,0.95)
)

rownames(LIP) <- c("AT","BE","CZ","EE","FI","FR","DE","GR","HU",
                   "IE","IT","NL","PL","PT","RO","ES","SE","UK")

superSubset(LIP, outcome = "S", incl.cut = 0.9)                   

testTESA(LIP, outcome = "S", incl.cut1 = 0.75)                   
}

