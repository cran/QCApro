\name{d.napoleon}

\docType{data}

\alias{d.napoleon}

\title{Reception of the Code Napoleon in Germany}

\description{
This dataset is from \href{http://dx.doi.org/10.1111/j.1748-121X.2009.00150.x}{Arvind and Stirton (2010)}, who analyze the reception of the
Code Napoleon in Germany using fsQCA.
}

\usage{data(d.napoleon)}

\format{
This data frame contains 14 rows (cases) and the following 8 columns (factors):

\tabular{rllrll}{
[ , 1]  \tab \bold{D} \tab exogenous factor\tab: \tab legal system              \tab ("1" heterogenous, "0" homogenous)\cr
[ , 2]  \tab \bold{C} \tab exogenous factor\tab: \tab territory                 \tab ("1" ruled by France, "0" ruled by enemy)\cr
[ , 3]  \tab \bold{I} \tab exogenous factor\tab: \tab state institutions        \tab ("1" strong, "0" none)\cr
[ , 4]  \tab \bold{F} \tab exogenous factor\tab: \tab economy                   \tab ("1" seigneural-feudal, "0" proto-industrial)\cr
[ , 5]  \tab \bold{L} \tab exogenous factor\tab: \tab ideology of state ruler   \tab ("1" liberal, "0" conservative)\cr
[ , 6]  \tab \bold{N} \tab exogenous factor\tab: \tab nativist tendencies       \tab ("1" yes, "0" no)\cr
[ , 7]  \tab \bold{A} \tab exogenous factor\tab: \tab sentiments towards France \tab ("1" very negative, "0" very positive)\cr
[ , 8]  \tab \bold{O} \tab endogenous factor\tab:   \tab adoption of Code Napoleon \tab ("1" yes, "0" no)
}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Arvind, Thiruvallore T., and Lindsay Stirton. 2010. \dQuote{Explaining the 
Reception of the Code Napoleon in Germany: A Fuzzy-Set Qualitative Comparative 
Analysis.} \emph{Legal Studies} \strong{30} (1):1-29.
DOI: \href{http://dx.doi.org/10.1111/j.1748-121X.2009.00150.x}{10.1111/j.1748-121X.2009.00150.x}.
}

\keyword{datasets}
