\name{superSubset}
\alias{superSubset}
\alias{is.sS}

\title{Find Superset and Subset Relations}

\description{
This helper function finds all combinations of conditions among all possible 
combinations that optimize the fulfilment of the specified criteria for a superset (necessity) or subset (sufficiency) relation to the outcome.
}

\usage{
superSubset(data, outcome = "", neg.out = FALSE, exo.facs = c(""), 
            relation = "nec", incl.cut = 1, cov.cut = 0, 
            use.tilde = FALSE, use.letters = FALSE, ...)

is.sS(x)
}

\arguments{
  \item{data}{A dataset of bivalent or multivalent crisp-set factor or bivalent 
              fuzzy-set variables.}
  \item{outcome}{The name of the outcome.}
  \item{neg.out}{Logical, use negation of \code{outcome}.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{relation}{The relation to \code{outcome}, either \code{"nec"},
        \code{"suf"}, \code{"necsuf"} or \code{"sufnec"}.}
  \item{incl.cut}{The minimal inclusion score of the relation.}
  \item{cov.cut}{The minimal coverage score of the relation.}
  \item{use.tilde}{Logical, use "~" for negation with bivalent variables.}
  \item{use.letters}{Logical, use simple letters instead of original factor
        names.}
  \item{...}{Other arguments for backward compatibility.}
  \item{x}{An object of class "sS".}
}

\details{
This helper function to the \code{testTESA} function returns a list of those of the 
\eqn{\prod_{j = 1}^{k}{(p_{j} + 1) - 1}}{\prod_{j = 1}^{k}{(p_{j} + 1) - 1}} 
potential value combinations, where \eqn{p_{j}}{p_{j}} is the number of values 
for exogenous variable \eqn{j}{j} and \eqn{k}{k} is the number of exogenous 
variables, that define minimal condition sets for the specified inclusion 
(consistency) and coverage score cut-offs with respect to an outcome.

If \code{relation = "nec"} (default), the function finds (combinations of) 
conditions that are supersets of (necessary for) the outcome. It starts with an 
initiation set, which is comprised of all \eqn{\sum_{j = 1}^{k}{p_{j}}} simple 
condition sets. This set is expanded by incrementally forming set-theoretic 
intersections of a higher order as long as \code{incl.cut} and \code{cov.cut} 
are still met (the former always takes precedence over the latter). If suitable
conjunctions exist, they will be returned, together with all their lower-order
conjuncts.

If none of the simple conditions or their negations in the initiation set passes 
\code{incl.cut}, disjunctions instead of conjunctions are formed until 
\code{incl.cut} and \code{cov.cut} will have been met. Only the disjunctions thus 
found will be returned.

If \code{relation = "suf"}, the function finds (combinations of) 
conditions that are subsets of (sufficient for) the outcome. The initiation set 
is comprised of all \eqn{\prod_{j = 1}^{k}{p_{j}}}{\prod_{j = 1}^{k}{p_{j}}} 
intersections of order \eqn{k}{k}. This set is reduced by incrementally forming 
intersections of a lower order as long as \code{incl.cut} and \code{cov.cut} are 
still met. Only the intersections of the lowest order will be printed. For more 
details, see Thiem and Dusa (2013). For \code{relation = "necsuf"} and 
\code{relation = "sufnec"}, \code{incl.cut} will be applied to each relation and 
\code{cov.cut} has no effect.

The argument \code{outcome} specifies the outcome. Outcomes from multivalent 
variables require curly-bracket notation (\code{X{value}}).

The logical argument \code{neg.out} controls whether \code{outcome} is to be 
used or its negation. If \code{outcome} is from a multivalent crisp-set factor, 
\code{neg.out = TRUE} has the effect that the disjunctions of all remaining values 
becomes the new outcome.

The argument \code{exo.facs} specifies the names of the exogenous factors. 
If omitted, all factors in \code{data} are used except the factor of which 
\code{outcome} is a level.

The argument \code{use.tilde} only applies to bivalent factors. If factors are 
already named with single letters, the argument \code{use.letters} has no effect.
}

\value{An object of class 'sS', which is a list with the following two main components:\cr
  \item{incl.cov}{A data frame with the parameters of fit.}
  \item{coms}{A data frame with the combination membership scores.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}
    
\references{
Ragin, Charles C. 2009. \dQuote{Qualitative Comparative Analysis Using Fuzzy Sets (fsQCA).} In \emph{Configurational Comparative Methods: Qualitative Comparative Analysis (QCA) and Related Techniques}, ed. B. Rihoux and C. C. Ragin. London: Sage Publications, pp. 87-121.

Schneider, Carsten Q., and Claudius Wagemann. 2013. \dQuote{Doing Justice to 
Logical Remainders in QCA: Moving Beyond the Standard Analysis.} \emph{Political 
Research Quarterly} \strong{66} (1):211-20. DOI: \href{http://dx.doi.org/10.1177/1065912912468269}{10.1177/1065912912468269}.

Thiem, Alrik. 2015. \emph{Standards of Good Practice and the Methodology of Necessary Conditions in Qualitative Comparative Analysis: A Critical View on Schneider and Wagemann's Theory-Guided/Enhanced Standard Analysis}. COMPASSS WP Series 2015-83. URL: \href{http://www.compasss.org/wpseries/Thiem2015.pdf}{http://www.compasss.org/wpseries/Thiem2015.pdf}.   
}

\seealso{\code{\link{testTESA}}}

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

# with mv-data from Hartmann and Kemmerzell (2010)
#-------------------------------------------------
data(d.partybans)
head(d.partybans)

HK <- superSubset(d.partybans, outcome = "PB", 
  exo.facs = c("C", "F", "T", "V"), incl.cut = 0.75)
HK

# combination membership scores for all cases (only first four 
# combinations and first ten lines displayed)
HK$coms[1:10, 1:4 , drop = FALSE]
}

\keyword{functions}
