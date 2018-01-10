\name{ambiguity}

\alias{ambiguity}

\title{Analyze the Degree of Ambiguity across Research Design Variations}

\description{
This evaluation function computes the degree of ambiguity across variations of a
reference research design. It has initially been programmed for \href{http://dx.doi.org/10.1177/0049124115610351}{Baumgartner and Thiem (2017)}.
}

\usage{
ambiguity(data, outcome = c(""), neg.out = c(FALSE), exo.facs = c(""),
          tuples = c(), incl.cut1 = c(1), incl.cut0 = c(1), sol.type = c("ps"), 
          row.dom = c(FALSE), min.dis = c(FALSE))
}

\arguments{
  \item{data}{A set of configurational data as processable by the 
        \code{\link{eQMC}} function.}
  \item{outcome}{A character vector of outcomes.}
  \item{neg.out}{A logical vector specifying whether to negate outcomes.}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{tuples}{A numeric vector of tuples of exogenous factors to be created 
        from \code{exof.cols}.}
  \item{incl.cut1}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{incl.cut0}{The maximum sufficiency inclusion score for an output function 
        value of "0".}
  \item{sol.type}{A character vector specifying the solution types to be generated.}
  \item{row.dom}{A logical vector imposing row dominance as a constraint on the solution to
        eliminate dominated inessential prime implicants.}
  \item{min.dis}{A logical vector imposing minimal disjunctivity as a constraint 
        on the solution to eliminate models with more prime implicants than the 
        model(s) with the fewest prime implicants.}
}

\details{
This evaluation function computes the degree of ambiguity across variations of a
reference design by recording the number of models for each design solution. It has initially been programmed for \href{http://dx.doi.org/10.1177/0049124115610351}{Baumgartner and Thiem (2015)}.

The argument \code{data} requires a set of configurational data as processable by 
the \code{\link{eQMC}} function.  

The argument \code{outcome} is a character vector, specifying the outcome(s) to be analyzed, either in curly-bracket notation (e.g., \code{O{value}}) if the outcome is from a multivalent (or a bivalent) factor, or in upper-case notation if the outcome is from a bivalent factor (e.g., \code{O} as a short-cut for \code{O{1}}). Outcomes from multivalent crisp-set factors always require curly-bracket notation. Outcomes can be single levels of factors not simultaneously passed to \code{exo.facs}. At least one outcome has to be specified.

The argument \code{neg.out} requires a logical vector of length one or two, whose values, which must not be duplicated, specify whether to negate the outcomes determined by \code{outcome}. If an element in \code{outcome} is a level from a multivalent factor, \code{neg.out = TRUE} makes the disjunction of all remaining levels the outcome. Possible values for \code{neg.out} include \code{FALSE}, \code{TRUE}, \code{FALSE, TRUE} and \code{TRUE, FALSE}.

The argument \code{exo.facs} is a character vector with the names of the exogenous factors. If omitted, all factors in \code{data} are used except that/those of the outcome/s given in \code{outcome}. and \code{tuples} specifies a numeric vector of tuples of exogenous factors to be created from \code{exo.facs}.

Minterms with an inclusion score of at least \code{incl.cut1} are coded positive (\code{OUT = "1"}), minterms with an inclusion score below \code{incl.cut1} but with at least \code{incl.cut0} are coded as a contradiction (\code{OUT = "C"}), and minterms with an inclusion score below \code{incl.cut0} are coded negative (\code{OUT = "0"}). If \code{inc .cut0} is not explicitly changed, it is set equal to \code{incl.cut1}.

The argument \code{sol.type} requires a character vector specifying the solution types to be generated. For example, \code{c("ps", "cs")} means parsimonious and conservative solution type.

The argument \code{row.dom} requires a logical vector, and controls whether the principle of row dominance is imposed as a constraint on the solution. An inessential prime implicant \eqn{P} dominates another \eqn{Q} if all configurations covered by \eqn{Q} are also covered by \eqn{P}, but they are not interchangeable (cf. \href{http://dx.doi.org/10.1002/j.1538-7305.1956.tb03835.x}{McCluskey 1956, 1425}; McCluskey 1965, 164-152). If row dominance is operative, models that contain dominated prime implicants will not be returned.

The argument \code{min.dis} requires a logical vector, and controls whether the principle of minimal disjunctivity is imposed as a constraint on the solution (McCluskey 1965, 12 -126). If minimal disjunctivity is operative, models that contain more than the number of prime implicants of the model(s) with the fewest prime implicants will not be returned.
}

\value{A list with the following two main components:\cr
  \item{tuples}{A list of all tuples of exogenous factors of the respective size 
        taken from all factors given in \code{exo.facs}.}
  \item{n.models}{A list of matrices giving the number of models in each solution
        for each design. The coding of labels has the following structure: 
        \code{O.1234}, where \code{O} is the outcome, \code{1} specifies the value 
        at the respective index of the argument \code{neg.out}, \code{2} the value 
        at the respective common index of the arguments \code{incl.cut1} and 
        \code{incl.cut0}, \code{3} the value at the respective index of the argument
        \code{sol.type}, and \code{4} the value at the respective common index of 
        the arguments \code{row.dom} and \code{min.dis}.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors (alphabetical)}{\tabular{ll}{
Thiem, Alrik        \tab: development, documentation, programming, testing
}}

\references{
Baumgartner, Michael, and Alrik Thiem. 2017. \dQuote{Model Ambiguities in Configurational Comparative Research.} \emph{Sociological Methods & Research} \strong{46} (4):954-87. DOI: \href{http://dx.doi.org/10.1177/0049124115610351}{10.1177/0049124115610351}.

McCluskey, Edward J. 1956. \dQuote{Minimization of Boolean Functions.} 
\emph{Bell Systems Technical Journal} \strong{35} (6):1417-44. DOI: \href{http://dx.doi.org/10.1002/j.1538-7305.1956.tb03835.x}{10.1002/j.1538-7305.1956.tb03835.x}.

McCluskey, Edward J. 1965. \emph{Introduction to the Theory of Switching Circuits}.
Princeton: Princeton University Press.
}

\seealso{\code{\link{eQMC}}}

\examples{
\dontrun{
# load dataset
data(d.tumorscreen)

# designs: outcomes HPF and LPF; all 3 to 5-tuples of exogenous factors
designs <- ambiguity(d.tumorscreen, outcome = c("HPF", "LPF"), 
  neg.out = c(FALSE, TRUE), tuples = 3:5)

# share of solutions with ambiguities
mapply(function (x) round(colSums((x > 1)) / nrow(x), 2), designs$n.models)
}
}

\keyword{functions}
