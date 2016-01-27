\name{d.represent}

\docType{data}

\alias{d.represent}

\title{Women's Representation in Western Democratic Parliaments}

\description{
This dataset is from \href{http://dx.doi.org/10.1111/j.1467-9248.2010.00833.x}{Krook (2010)}, who analyzes the determinants of high women's 
representation in Western-democratic parliaments using csQCA.
}

\usage{data(d.represent)}

\format{
This data frame contains 22 rows (cases) and the following 6 columns (factors):

\tabular{rlll}{
[ , 1] \tab \bold{ES}  \tab exogenous factor: \tab PR electoral system                              \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                \cr
[ , 2] \tab \bold{QU}  \tab exogenous factor: \tab quota for women                                  \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                \cr
[ , 3] \tab \bold{WS}  \tab exogenous factor: \tab social-democratic welfare system                 \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                \cr
[ , 4] \tab \bold{WM}  \tab exogenous factor: \tab autonomous women's movement                      \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                \cr
[ , 5] \tab \bold{LP}  \tab exogenous factor: \tab seats held by left-libertarian parties           \cr
       \tab            \tab            \tab ("1" at least 7 percent, "0" less than 7 percent)\cr
[ , 6] \tab \bold{WNP} \tab endogenous factor:   \tab women in single/lower house of parliament        \cr
       \tab            \tab            \tab ("1" at least 30 percent, "0" less than 30 percent)
}
}

\note{
The row names are the official International Organization for Standardization (ISO) country code elements as 
specified in ISO 3166-1-alpha-2.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Krook, Mona Lena. 2010. \dQuote{Women's Representation in Parliament: A Qualitative 
Comparative Analysis.} \emph{Political Studies} \strong{58} (5):886-908.
DOI: \href{http://dx.doi.org/10.1111/j.1467-9248.2010.00833.x}{10.1111/j.1467-9248.2010.00833.x}.
}


\keyword{datasets}