\name{d.tumorscreen}

\docType{data}

\alias{d.tumorscreen}

\title{Comparing Universal Lynch Syndrome Tumor-Screening Programs}

\description{
This dataset is from \href{http://dx.doi.org/10.1177/1078087403262861}{Cragun \emph{et al}. (2014)}, who analyze the association between different universal tumor screening procedures and certain levels of patient follow-through with germ-line testing for Lynch Syndrome after a screen-positive result using csQCA.
}

\usage{data(d.tumorscreen)}

\format{
This data frame contains 15 rows (cases) and the following 8 columns (factors):

\tabular{rllrll}{
[ , 1]  \tab \bold{HPF} \tab endogenous factor\tab: \tab high patient follow-through \tab ("1" yes, "0" no)\cr
[ , 2]  \tab \bold{LPF} \tab endogenous factor\tab: \tab low patient follow-through      \tab ("1" yes, "0" no)\cr
[ , 3]  \tab \bold{CA}  \tab exogenous factor\tab: \tab challenges to adoption at least as high as facilitators       \tab ("1" yes, "0" no)\cr
[ , 4]  \tab \bold{AR} \tab exogenous factor\tab: \tab automatic reflex test of screen-positive tumors \tab ("1" yes, "0" no)\cr
[ , 5]  \tab \bold{RR}  \tab exogenous factor\tab: \tab genetic counselor receives positive screen results                 \tab ("1" yes, "0" no)\cr
[ , 6]  \tab \bold{DR}  \tab exogenous factor\tab:   \tab genetic counselor discloses screening result to patient               \tab ("1" yes, "0" no)\cr
[ , 7]  \tab \bold{DC}  \tab exogenous factor\tab:   \tab difficulty in contacting patients               \tab ("1" yes, "0" no)\cr
[ , 8]  \tab \bold{PR}  \tab exogenous factor\tab:   \tab need for physician referral is a barrier               \tab ("1" yes, "0" no)
}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Cragun, Deborah, Rita D. DeBate, Susan T. Vadaparampil, Julie Baldwin, Heather Hampel, and Tuya Pal. 2014. \dQuote{Comparing Universal Lynch Syndrome Tumor-Screening Programs to Evaluate Associations between Implementation Strategies and Patient Follow-Through.} \emph{Genetics in Medicine} \strong{16} (10):773-82.
DOI: \href{http://dx.doi.org/10.1038/gim.2014.31}{10.1038/gim.2014.31}.
}

\keyword{datasets}
