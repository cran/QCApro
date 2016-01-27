\name{d.transport}

\docType{data}

\alias{d.transport}

\title{Determinants of High Transport Project Acceptance}

\description{
This dataset is from \href{http://dx.doi.org/10.1177/1098214011411574}{Sager and Andereggen (2012)}, who analyze the determinants 
of high transport project acceptance in Switzerland using mvQCA. 
}

\usage{data(d.transport)}

\format{
This data frame contains 21 rows (cases) and the following 10 columns (factors):

\tabular{rlll}{
[ , 1]  \tab \bold{FED} \tab exogenous factor: \tab federal level                             \cr
        \tab            \tab            \tab ("2" federal, "1" cantonal, "0" municipal)\cr
[ , 2]  \tab \bold{FIN} \tab exogenous factor: \tab financial situation                       \cr
        \tab            \tab            \tab ("1" positive, "0" negative)              \cr
[ , 3]  \tab \bold{URB} \tab exogenous factor: \tab sociostructural project location          \cr
        \tab            \tab            \tab ("1" urban, "0" rural)                    \cr
[ , 4]  \tab \bold{GER} \tab exogenous factor: \tab cultural project location                 \cr
        \tab            \tab            \tab ("1" German-speaking, "0" French-speaking)\cr
[ , 5]  \tab \bold{HIS} \tab exogenous factor: \tab prior history                             \cr
        \tab            \tab            \tab ("1" yes, "0" no)                         \cr
[ , 6]  \tab \bold{COO} \tab exogenous factor: \tab planning coordination                     \cr
        \tab            \tab            \tab ("1" strong, "0" not strong)              \cr
[ , 7]  \tab \bold{PRO} \tab exogenous factor: \tab administrative professionalization        \cr
        \tab            \tab            \tab ("1" high, "0" not high)                  \cr
[ , 8]  \tab \bold{DIS} \tab exogenous factor: \tab administration's discretion               \cr
        \tab            \tab            \tab ("1" broad, "0" not broad)                \cr
[ , 9]  \tab \bold{EXP} \tab exogenous factor: \tab influence of external experts             \cr
        \tab            \tab            \tab ("1" great, "0" not great)                \cr
[ , 10] \tab \bold{ACC} \tab endogenous factor:   \tab project acceptance                        \cr
        \tab            \tab            \tab ("1" high, "0" not high)
}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Sager, Fritz, and Celine Andereggen. 2012. \dQuote{Dealing With Complex Causality 
in Realist Synthesis: The Promise of Qualitative Comparative Analysis.} 
\emph{American Journal of Evaluation} \strong{33} (1):60-78.
DOI: \href{http://dx.doi.org/10.1177/1098214011411574}{10.1177/1098214011411574}.
}

\keyword{datasets}
