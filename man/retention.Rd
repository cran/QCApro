\name{retention}
\alias{retention}

\title{Compute Retention Probabilities of QCA Baseline Solutions}

\description{
This evaluation function computes retention probabilities of QCA baseline solutions. It has 
been programmed for \href{http://dx.doi.org/10.1093/pan/mpv028}{Thiem, Spoehel, and Dusa (2016)}.
}

\usage{
retention(data, outcome = "", exo.facs = c(""), type = "corruption",
          assump = "DPA", n.cut = 1, incl.cut = 1, p.pert = 0.5,
          n.pert = 1)
}

\arguments{
  \item{data}{A dataset of bivalent crisp-set factors.}
  \item{outcome}{The name of the outcome.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{type}{Induce errors on the endogenous factor or delete cases.}
  \item{assump}{Assume dependent or independent perturbations.}
  \item{n.cut}{The minimum number of cases for a minterm not to be 
        considered as a remainder.}
  \item{incl.cut}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{p.pert}{The probability of perturbation under independence.}
  \item{n.pert}{The number of perturbations under dependence.}
}

\details{
This function computes exact retention probabilities of QCA baseline solutions for saturated truth tables and truth tables with a two-difference restriction (every remainder differs on at least two positions from every positive minterm). 

The argument \code{data} requires a suitable dataset. Suitable datasets have the following structure: values of "0" and "1" for bivalent crisp-set factors. 

The argument \code{exo.facs} specifies the exogenous factors. If omitted, all 
factors in \code{data} are used except that of the \code{outcome}.

The argument \code{type} specifies whether errors are to be induced in the endogenous factor ("1" is recoded to "0"; "0" is recoded to "1") of cases or whether entire cases are to be deleted from the data.

The argument \code{assump} specifies whether the perturbations detailed in \code{type} occur independently of each other or whether they are dependent on each other. Note that the assumption of dependence increases the consumption of computational resources significantly.

Minterms that contain fewer than \code{n.cut} cases with membership scores above 0.5 are coded as remainders (\code{OUT = "?"}). If the number of such cases is at least \code{n.cut}, minterms with an inclusion score of at least \code{incl.cut} are coded positive (\code{OUT = "1"}), and minterms with an inclusion score below \code{incl.cut} are coded negative (\code{OUT = "0"}). The possibility to specify contradictions using a second inclusion cut-off as in the \code{\link{truthTable}} function does not exist.

The argument \code{p.pert} specifies the probability of perturbation for \code{type = "independent"}. For example, if \code{p.pert = 1}, each case is guaranteed to have measurement error on the endogenous factor.

The argument \code{n.pert} specifies the number of perturbations for \code{type = "dependent"}. This must be an integer between zero (no case suffers from measurement error in the endogenous factor or no case gets deleted) and the total number of cases in \code{data} (all cases suffer from measurement error in the endogenous factor or all cases get deleted.)
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: programming, testing\cr
Spoehel, Reto\tab: development\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Hug, Simon. 2013. \dQuote{Qualitative Comparative Analysis: How Inductive Use and 
Measurement Error lead to Problematic Inference.} \emph{Political Analysis} \strong{21} (2):252-65.
DOI: \href{http://dx.doi.org/10.1093/pan/mps061}{10.1093/pan/mps061}.

Thiem, Alrik, Reto Spoehel, and Adrian Dusa. 2016. \dQuote{Enhancing Sensitivity Diagnostics for Qualitative Comparative Analysis: A Combinatorial Approach.} \emph{Political Analysis} \strong{24} (1):104-20. DOI: \href{http://dx.doi.org/10.1093/pan/mpv028}{10.1093/pan/mpv028}.
}

\seealso{\code{\link{truthTable}}}

\examples{
# replicate results from Hug (2013) for 2 deleted cases
#------------------------------------------------------
dat <- data.frame(matrix(c(
   rep(1,25),rep(0,20),rep(c(0,0,1,0,0),3),
   0,0,0,1,0,0,1,0,0,0,0,rep(1,7),0,1),
   nrow = 16, byrow = TRUE, dimnames = list(c(
    "AT","DK","FI","NO","SE","AU","CA","FR",
    "US","DE","NL","CH","JP","NZ","IE","BE"),
   c("P", "U", "C", "S", "W"))
))

retention(dat, outcome = "W", type = "deletion", assump = "dependent", n.pert = 2)
}

\keyword{functions}
