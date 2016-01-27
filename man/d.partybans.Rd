\name{d.partybans}

\docType{data}

\alias{d.partybans}

\title{Party Bans Sub-Saharan Africa}

\description{
This dataset is from \href{http://dx.doi.org/10.1080/13510347.2010.491189}{Hartmann and Kemmerzell (2010)}, who analyze the determinants of
the introduction of party ban provisions and their actual implementation in sub-Saharan Africa using mvQCA.
}

\usage{data(d.partybans)}

\format{
This data frame contains 48 rows (cases) and the following 7 columns (factors):

\tabular{rlll}{
[ , 1] \tab \bold{C}   \tab exogenous factor: \tab colonial tradition                                   \cr
       \tab            \tab            \tab ("2" British, "1" French, "0" other)                 \cr
[ , 2] \tab \bold{F}   \tab exogenous factor: \tab former regime type competition                       \cr
       \tab            \tab            \tab ("2" no, "1" limited, "0" multi-party)               \cr
[ , 3] \tab \bold{T}   \tab exogenous factor: \tab mode of transition                                   \cr
       \tab            \tab            \tab ("2" managed, "1" pacted, "0" democracy before 1990) \cr
[ , 4] \tab \bold{R}   \tab exogenous factor: \tab regime type                                          \cr
       \tab            \tab            \tab ("2" authoritarian, "1" liberalizing, "0" democratic)\cr
[ , 5] \tab \bold{V}   \tab exogenous factor: \tab ethnic violence                                      \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                    \cr
[ , 6] \tab \bold{PB}  \tab endogenous factor:   \tab party ban provisions introduced                      \cr
       \tab            \tab            \tab ("1" yes, "0" no)                                    \cr
[ , 7] \tab \bold{PBI} \tab endogenous factor:   \tab party bans implemented                               \cr
       \tab            \tab            \tab ("1" yes, "0" no)
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
Hartmann, Christof, and Joerg Kemmerzell. 2010. \dQuote{Understanding Variations 
in Party Bans in Africa.} \emph{Democratization} \strong{17} (4):642-65.
DOI: \href{http://dx.doi.org/10.1080/13510347.2010.491189}{10.1080/13510347.2010.491189}.
}

\keyword{datasets}