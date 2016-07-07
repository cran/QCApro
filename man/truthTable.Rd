\name{truthTable}
\alias{truthTable}
\alias{is.tt}

\title{Create a Truth Table}

\description{
This function creates truth tables from configurational data.
}

\usage{
truthTable(data, outcome = "", neg.out = FALSE, exo.facs = c(""), 
           n.cut = 1, incl.cut1 = 1, incl.cut0 = 1, complete = FALSE, 
           show.cases = FALSE, sort.by = c(""), decreasing = TRUE, 
           inf.test = c(""), use.letters = FALSE, ...)

is.tt(x)
}

\arguments{
  \item{data}{A set of configurational data (of class 'matrix' or 'data.frame').}
  \item{outcome}{The name of the outcome.}
  \item{neg.out}{Logical, use the negation of \code{outcome}.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{n.cut}{The minimum number of cases with set membership score above 0.5 
        for an output function value of "0", "1" or "C"; an integer between 1 and
        the maximum number of cases for all non-remainder minterms.}
  \item{incl.cut1}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{incl.cut0}{The maximum sufficiency inclusion score for an output function 
        value of "0".}
  \item{complete}{Logical, print the complete truth table.}
  \item{show.cases}{Logical, print case names (do not use this option with many 
        cases and/or long case names).}
  \item{sort.by}{Sort the truth table by inclusion scores and/or number of cases.}
  \item{decreasing}{Sort in decreasing or increasing order of value(s) passed to 
        \code{sort.by}.}
  \item{inf.test}{A vector of length two specifying the inference-statistical 
        test to be performed (currently only \code{"binom"}) and the critical
        significance level.}
  \item{use.letters}{Logical, use single letters (in alphabetical order) instead 
        of original variable names.}
  \item{...}{Other arguments.}
  \item{x}{An object of class 'tt'.}
}

\details{
The argument \code{data} can be a truth table object (an object of class 'tt' 
returned by the \code{truthTable} function) or a suitable data set. Suitable data
sets have the following structure: values of 0 and 1 for bivalent crisp-set factors, 
values between 0 and 1 for bivalent fuzzy-set factors, and values beginning 
with 0 at increments of 1 for multivalent crisp-set factors. The placeholders 
"-" and "dc" indicate "don't cares" in auxiliary factors that specify temporal 
order between other substantive factors in tQCA. These values lead to the 
exclusion of the auxiliary factor from the computation of parameters of fit.

The argument \code{outcome} specifies the outcome to be analyzed, either in
curly-bracket notation (e.g., \code{O{value}}) if the outcome is from a multivalent
(or a bivalent) factor, or in upper-case notation if the outcome is from a bivalent
factor (e.g., \code{O} as a short-cut for \code{O{1}}). Outcomes from multivalent crisp-set factors always require curly-bracket notation. Outcomes must be single 
levels of factors not simultaneously passed to \code{exo.facs}.

The logical argument \code{neg.out} controls whether \code{outcome} is to be 
analyzed or its negation. If \code{outcome} is a level from a multivalent factor, 
\code{neg.out = TRUE} causes the disjunction of all remaining levels to become 
the outcome to be analyzed.

The argument \code{exo.facs} specifies the exogenous factors. If omitted, all 
factors in \code{data} are used except that of the \code{outcome}. Please note 
that computation times may increase significantly beyond 17 exogenous factors, 
and that the computation of a solution may not be possible at all depending on 
end-user machine constraints.

Minterms that contain fewer than \code{n.cut} cases with membership scores above 
0.5 are coded as remainders (\code{OUT = "?"}). If the number of such cases is at 
least \code{n.cut}, minterms with an inclusion score of at least \code{incl.cut1} 
are coded positive (\code{OUT = "1"}), minterms with an inclusion score below 
\code{incl.cut1} but with at least \code{incl.cut0} are coded as a contradiction 
(\code{OUT = "C"}), and minterms with an inclusion score below \code{incl.cut0} 
are coded negative (\code{OUT = "0"}). If \code{incl.cut0} is not explicitly 
changed, it is set equal to \code{incl.cut1}.

The logical argument \code{show.cases} controls whether case names are displayed
next to their corresponding minterm (do not use this option with many cases 
and/or long case names).

The \code{sort.by} argument orders all minterms by their inclusion scores 
(\emph{incl}) or the number of cases with membership above 0.5 they contain 
(\emph{n}) or both, in either order.

If the exogenous factors are already named with single letters, the argument 
\code{use.letters} will have no effect when set to \code{TRUE}. Otherwise, 
upper-case letters will replace original factor names in alphabetical order.

The argument \code{inf.test} provides functionality for basing output function
value codings on inference-statistical tests. Currently, only an exact binomial
test (\code{"binom"}) is available, which requires the data to contain only
bivalent or multivalent crisp-set factors. The argument requires a vector of
length two, comprising the test and a critical significance level. If the
empirical inclusion score of a minterm is not significantly lower than 
\code{incl.cut1}, it will be coded positive (\code{OUT = "1"}). If it is 
significantly lower than \code{incl.cut1} yet still significantly higher than 
\code{incl.cut0}, it will be coded as a contradiction (\code{OUT = "C"}). If it
is not significantly higher than \code{incl.cut0}, it will be coded negative
(\code{OUT = "0"}).
}

\value{An object of class 'tt', which is a list with the following six main 
       components:\cr
  \item{tt}{The truth table.}
  \item{indexes}{The minterm line numbers.}
  \item{noflevels}{A vector with the number of levels of the exogenous factors.}
  \item{initial.data}{The initial data.}
  \item{recoded.data}{Recoded data (if crisp, same as \code{initial.data}; if 
                      fuzzy, dichotomized version of \code{initial.data}).}
  \item{cases}{The cases with membership above 0.5 in a minterm.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming\cr
Thiem, Alrik\tab: development, documentation, programming, testing
}}

\references{
Emmenegger, Patrick. 2011. \dQuote{Job Security Regulations in Western Democracies: 
A Fuzzy Set Analysis.} \emph{European Journal of Political Research} \strong{50} (3):336-64.
DOI: \href{http://dx.doi.org/10.1111/j.1475-6765.2010.01933.x}{10.1111/j.1475-6765.2010.01933.x}.

Hartmann, Christof, and Joerg Kemmerzell. 2010. \dQuote{Understanding Variations 
in Party Bans in Africa.} \emph{Democratization} \strong{17} (4):642-65.
DOI: \href{http://dx.doi.org/10.1080/13510347.2010.491189}{10.1080/13510347.2010.491189}.

Krook, Mona Lena. 2010. \dQuote{Women's Representation in Parliament: A Qualitative 
Comparative Analysis.} \emph{Political Studies} \strong{58} (5):886-908.
DOI: \href{http://dx.doi.org/10.1111/j.1467-9248.2010.00833.x}{10.1111/j.1467-9248.2010.00833.x}.

Ragin, Charles C., and Sarah Ilene Strand. 2008. \dQuote{Using Qualitative 
Comparative Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).} 
\emph{Sociological Methods & Research} \strong{36} (4):431-41.
DOI: \href{http://dx.doi.org/10.1177/0049124107313903}{10.1177/0049124107313903}.
}

\examples{
# csQCA using Krook (2010)
#-------------------------
data(d.represent)
head(d.represent)

# print truth table; if all factors except that of the outcome in 
# the data should be included as exogenous factors, then these need 
# not be specified separately
truthTable(d.represent, outcome = "WNP")

# print complete truth table, show cases, and first sort by 
# inclusion scores, then by number of cases
truthTable(d.represent, outcome = "WNP", complete = TRUE, 
  show.cases = TRUE, sort.by = c("incl", "n"))

# code minterms with a single case as remainders (note: use of 
# 'n.cut' should be well justified)
KRO.tt <- truthTable(d.represent, outcome = "WNP", n.cut = 2, 
  show.cases = TRUE)
KRO.tt

# print cases that were assigned to remainders based on argument 'n.cut'
KRO.tt$excluded

# fsQCA using Emmenegger (2011)
#------------------------------
data(d.jobsecurity)
head(d.jobsecurity)

# code non-remainder minterms with inclusion scores between 0.4 
# and 0.8 as contradictions (note: these are not 'contradictions' 
# in the logical sense of the word but minterms that can neither 
# be coded as sufficient nor as insufficient for the outcome)
truthTable(d.jobsecurity, outcome = "JSR", incl.cut1 = 0.8, incl.cut0 = 0.4)

# truth table based on the negated outcome
truthTable(d.jobsecurity, outcome = "JSR", neg.out = TRUE, incl.cut1 = 0.8, 
  incl.cut0 = 0.4)

# mvQCA using Hartmann and Kemmerzell (2010)
#-------------------------------------------
data(d.partybans)
head(d.partybans)

# code non-remainder minterms with inclusion scores below 1 
# but above 0.4 as contradictions 
HK.tt <- truthTable(d.partybans, outcome = "PB", 
  exo.facs = c("C","F","T","V"), incl.cut0 = 0.4)
HK.tt

# list the number of levels for the exogenous factors  
HK.tt$noflevels

# which minterms have more than 2 cases?
HK.tt$tt[which(HK.tt$tt$n > 2), ]

# code output function values in truth table based on 
# exact binomial test
truthTable(d.partybans, outcome = "PB", exo.facs = c("C","F","T"), 
  incl.cut1 = 0.9, incl.cut0 = 0.4, show.cases = TRUE, 
  inf.test = c("binom", 0.1))

# tQCA using Ragin and Strand (2008)
#-----------------------------------
data(d.graduate)
head(d.graduate)

# tQCA truth table with "don't care" values
truthTable(d.graduate, outcome = "REC")

}

\keyword{functions}
