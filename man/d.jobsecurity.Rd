\name{d.jobsecurity}

\docType{data}

\alias{d.jobsecurity}

\title{Job Security Regulations in Western Democracies}

\description{
This dataset is from \href{http://dx.doi.org/10.1111/j.1475-6765.2010.01933.x}{Emmenegger (2011)}, who analyzes the determinants of 
high job security regulations in Western democracies using fsQCA. 
}

\usage{data(d.jobsecurity)}

\format{
This data frame contains 19 rows (cases) and the following 7 columns (factors):

\tabular{rllrll}{
[ , 1] \tab \bold{S}   \tab exogenous factor\tab: \tab level of statism                  \tab ("1" high, "0" not high)\cr
[ , 2] \tab \bold{C}   \tab exogenous factor\tab: \tab level of non-market coordination  \tab ("1" high, "0" not high)\cr
[ , 3] \tab \bold{L}   \tab exogenous factor\tab: \tab level of labour movement strength \tab ("1" high, "0" not high)\cr
[ , 4] \tab \bold{R}   \tab exogenous factor\tab: \tab level of Catholicism              \tab ("1" high, "0" not high)\cr
[ , 5] \tab \bold{P}   \tab exogenous factor\tab: \tab level of religious party strength \tab ("1" high, "0" not high)\cr
[ , 6] \tab \bold{V}   \tab exogenous factor\tab: \tab institutional veto points         \tab ("1" many, "0" not many)\cr
[ , 7] \tab \bold{JSR} \tab endogenous factor\tab:   \tab level of job security regulations \tab ("1" high, "0" not high)
}
}

\note{
The row names are the official International Organization for Standardization 
(ISO) country code elements as specified in ISO 3166-1-alpha-2.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Emmenegger, Patrick. 2011. \dQuote{Job Security Regulations in Western Democracies: 
A Fuzzy Set Analysis.} \emph{European Journal of Political Research} \strong{50} (3):336-64.
DOI: \href{http://dx.doi.org/10.1111/j.1475-6765.2010.01933.x}{10.1111/j.1475-6765.2010.01933.x}.
}

\keyword{datasets}
