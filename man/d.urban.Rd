\name{d.urban}

\docType{data}

\alias{d.urban}

\title{Explaining Urban Regimes}

\description{
This dataset is from \href{http://dx.doi.org/10.1177/1078087403262861}{Kilburn (2004)}, who analyzes the influence of city context
on urban regimes across 14 cities in the United States using csQCA.
}

\usage{data(d.urban)}

\format{
This data frame contains 14 rows (cases) and the following 6 columns (factors):

\tabular{rllrll}{
[ , 1]  \tab \bold{MLC} \tab exogenous factor\tab: \tab mobility of local capital \tab ("1" high, "0" not high)\cr
[ , 2]  \tab \bold{FRB} \tab exogenous factor\tab: \tab fiscal resource base      \tab ("1" large, "0" not large)\cr
[ , 3]  \tab \bold{CP}  \tab exogenous factor\tab: \tab civic participation       \tab ("1" high, "0" not high)\cr
[ , 4]  \tab \bold{WSR} \tab exogenous factor\tab: \tab ward-style representation \tab ("1" high, "0" not high)\cr
[ , 5]  \tab \bold{CS}  \tab exogenous factor\tab: \tab city size                 \tab ("1" large, "0" not large)\cr
[ , 6]  \tab \bold{RT}  \tab endogenous factor\tab:   \tab regime type               \tab ("1" progressive, "0" developmental/caretaker)
}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Kilburn, H. Whitt. 2004. \dQuote{Explaining U.S. Urban Regimes.}
\emph{Urban Affairs Review} \strong{39} (5):633-51.
DOI: \href{http://dx.doi.org/10.1177/1078087403262861}{10.1177/1078087403262861}.
}

\keyword{datasets}
