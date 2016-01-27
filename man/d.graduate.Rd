\name{d.graduate}

\docType{data}

\alias{d.graduate}

\title{Graduate Student Union Recognition at Research Universities}

\description{
This dataset is originally from \href{http://dx.doi.org/10.1177/0049124105277197}{Caren and Panofsky (2005)}, who analyze the 
determinants of unionization attempts by graduate student workers at research
universities using tQCA. Their study has been replicated and corrected by \href{http://dx.doi.org/10.1177/0049124107313903}{Ragin 
and Strand (2008)}.
}

\usage{data(d.graduate)}

\format{
This data frame contains 17 rows (cases) and the following 6 columns (factors):

\tabular{rllrlll}{
[ , 1] \tab \bold{P}   \tab exogenous factor\tab: \tab public university                \tab ("1" yes, "0" no)                \cr
[ , 2] \tab \bold{E}   \tab exogenous factor\tab: \tab support of elite allies          \tab ("1" yes, "0" no)                \cr
[ , 3] \tab \bold{A}   \tab exogenous factor\tab: \tab national union affiliation       \tab ("1" yes, "0" no)                \cr
[ , 4] \tab \bold{S}   \tab exogenous factor\tab: \tab a strike or strike threat        \tab ("1" yes, "0" no)                \cr
[ , 5] \tab \bold{EBA} \tab exogenous factor\tab: \tab \bold{E} present before \bold{A} \tab ("1" yes, "0" no, "-" don't care)\cr
[ , 6] \tab \bold{REC} \tab endogenous factor\tab:   \tab union recognition                \tab ("1" yes, "0" no)
}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{l}{
Thiem, Alrik: collection, documentation
}}

\source{
Caren, Neal, and Aaron Panofsky. 2005. \dQuote{TQCA: A Technique for Adding Temporality 
to Qualitative Comparative Analysis.} 
\emph{Sociological Methods & Research} \strong{34} (2):147-72.
DOI: \href{http://dx.doi.org/10.1177/0049124105277197}{10.1177/0049124105277197}.

Ragin, Charles C., and Sarah Ilene Strand. 2008. \dQuote{Using Qualitative Comparative 
Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).} 
\emph{Sociological Methods & Research} \strong{36} (4):431-41.
DOI: \href{http://dx.doi.org/10.1177/0049124107313903}{10.1177/0049124107313903}.
}

\keyword{datasets}