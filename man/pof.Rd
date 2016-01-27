\name{pof}

\alias{pof}
\alias{is.pof}

\title{Compute Set-Relational Parameters of Fit}

\description{
This function computes inclusion (consistency) and coverage scores.
}

\usage{
pof(setms, outcome, data, neg.out = FALSE, relation = "suf",
    inf.test = "", incl.cut1 = 0.75, incl.cut0 = 0.5, ...)

is.pof(x)
}

\arguments{
  \item{setms}{A data frame of set membership scores, or a matrix of implicants, 
        or a vector of implicant matrix line numbers.}
  \item{outcome}{The name of the outcome.}
  \item{data}{The working data set.}
  \item{neg.out}{Logical, use negation of \code{outcome}.}
  \item{relation}{The set relation to \code{outcome}, either \code{"nec"} or
        \code{"suf"}.}
  \item{inf.test}{The inference-statistical test to be performed (currently only
        \code{"binom"} for bivalent and multivalent crisp-set variables).}
  \item{incl.cut1}{The upper inclusion cut-off against which the
        empirical inclusion score is tested if \code{inf.test = "binom"}.}
  \item{incl.cut0}{The lower inclusion cut-off against which the
        empirical inclusion score is tested if \code{inf.test = "binom"}.}
  \item{...}{Other arguments (not used in this function).}
  \item{x}{An object of class "pof".}
}

\details{
The argument \code{setms} specifies a data frame of \emph{set} membership scores, 
where \emph{set} refers to any kind of set, including simple sets, combinations 
returned by the \code{superSubset} function (\code{coms}), prime implicants returned 
by the \code{eQMC} function (\code{pims}), or any other compound set.

The function also accepts a matrix of implicants with the level representation 
of created by the \code{mintermMatrix} function, or even a corresponding vector 
of implicant matrix line numbers.

The argument \code{outcome} specifies the outcome to be analyzed, either in
curly-bracket notation (e.g., \code{O{value}}) if the outcome is from a multivalent
(or a bivalent) factor, or in upper-case notation if the outcome is from a bivalent
factor (e.g., \code{O} as a short-cut for \code{O{1}}). Outcomes from multivalent crisp-set factors always require curly-bracket notation. Outcomes must be single 
levels of factors not simultaneously passed to \code{exo.facs}.

The logical argument \code{neg.out} controls whether \code{outcome} is to be 
analyzed or its negation. If \code{outcome} is a level from a multivalent factor, 
\code{neg.out = TRUE} causes the disjunction of all remaining levels to become 
the outcome to be analyzed.

The argument \code{inf.test} provides functionality for adjudicating between
rival hypotheses on the basis of inference-statistical tests. Currently, only an
exact binomial test (\code{"binom"}) is available, which requires the data to
contain only bivalent or multivalent crisp-set factors. Two one-tailed tests
are performed. The null hypothesis with respect to \code{incl.cut1} is that the
empirical inclusion score of each element in \code{setms} is not lower than the
upper critical inclusion cut-off provided in \code{incl.cut1}. The null hypothesis
with respect to \code{incl.cut0} is that the empirical inclusion score of each
element in \code{setms} is not higher than the lower critical inclusion cut-off
provided in \code{incl.cut0}.
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Emmenegger, Patrick. 2011. \dQuote{Job Security Regulations in Western Democracies: 
A Fuzzy Set Analysis.} \emph{European Journal of Political Research} \strong{50} 
(3):336-64. DOI: \href{http://dx.doi.org/10.1111/j.1475-6765.2010.01933.x}{10.1111/j.1475-6765.2010.01933.x}.

Hartmann, Christof, and Joerg Kemmerzell. 2010. \dQuote{Understanding Variations 
in Party Bans in Africa.} \emph{Democratization} \strong{17} (4):642-65.
DOI: \href{http://dx.doi.org/10.1080/13510347.2010.491189}{10.1080/13510347.2010.491189}.

Krook, Mona Lena. 2010. \dQuote{Women's Representation in Parliament: A Qualitative 
Comparative Analysis.} \emph{Political Studies} \strong{58} (5):886-908.
DOI: \href{http://dx.doi.org/10.1111/j.1467-9248.2010.00833.x}{10.1111/j.1467-9248.2010.00833.x}.
}

\seealso{\code{\link{eQMC}}}

\examples{
# csQCA using Krook (2010)
#-------------------------
data(d.represent)
head(d.represent)

# solution with details
KRO <- eQMC(d.represent, outcome = "WNP", incl.cut1 = 0.9, 
            details = TRUE)
KRO

# exact binomial tests of sufficiency inclusion
pof(KRO$pims, outcome = "WNP", d.represent, inf.test = c("binom", 0.1), 
    incl.cut1 = 0.75, incl.cut0 = 0.5)

# fsQCA using Emmenegger (2011)
#------------------------------
data(d.jobsecurity)
head(d.jobsecurity)

# solution with details
EMM.sol <- eQMC(d.jobsecurity, outcome = "JSR", incl.cut1 = 0.9, 
                details = TRUE)
EMM.sol

# are the model prime implicants also sufficient for the negation 
# of the outcome?
pof(EMM.sol$pims, outcome = "JSR", d.jobsecurity, neg.out = TRUE)

# are the negations of the model prime implicants also sufficient 
# for the outcome?
pof(1 - EMM.sol$pims, outcome = "JSR", d.jobsecurity)

# parameters of fit for matrix of implicants;
# "-1" is the placeholder for an eliminated variable;
# e.g.: R*p*V and S*c*L*P*v
#      "S"  "C"  "L"  "R"  "P"  "V"
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]   -1   -1   -1    1    0    1
#[2,]    1    0    1   -1    1    0
mat <- matrix(c(-1,-1,-1, 1, 0, 1, 
                 1, 0, 1,-1, 1, 0), nrow = 2, byrow = TRUE)

pof(mat, outcome = "JSR", d.jobsecurity)

# or even vectors of line numbers from the implicant matrix
pof(c(43, 57), "JSR", d.jobsecurity)


# mv-data from Hartmann and Kemmerzell (2010)
#-------------------------------------------
data(d.partybans)
head(d.partybans)

# parameters of fit for several mv-expressions
expr <- c("C{1}", "F{2}", "T{2}", "T{1}*V{0}")
dat <- data.frame(ifelse(d.partybans$C == 1, 1, 0),
                  ifelse(d.partybans$F == 2, 1, 0),
                  ifelse(d.partybans$T == 2, 1, 0), 
                  ifelse(d.partybans$T == 1 & d.partybans$V == 0, 1, 0))
colnames(dat) <- expr                        

pof(dat, outcome = "PB{1}", d.partybans)

# miscellaneous
#--------------
# parameters of fit for a data frame
x <- data.frame(A = c(1,1,1,0,1), B = c(1,1,1,0,1),
                C = c(0,1,0,0,1), D = c(0,0,1,0,1),
                O = c(1,1,1,0,1))

pof(x[, -5], outcome = "O", x)

# for a single column from that data frame
pof(x$A, x$O)

# for multiple columns from that data frame
pof(x[, 1:2], outcome = "O", x)

}

\keyword{functions}
