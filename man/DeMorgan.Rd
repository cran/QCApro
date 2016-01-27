\name{DeMorgan}

\alias{DeMorgan}
\alias{is.DeMorgan}

\title{Negate Boolean Expressions using De Morgan's Laws}

\description{
This function negates simple or complex Boolean expressions using the two 
De Morgan Laws.
}

\usage{
DeMorgan(expression, and.split = "", use.tilde = FALSE)

is.DeMorgan(x)
}

\arguments{
  \item{expression}{A string representing a Boolean expression or a solution 
                    object of class 'qca'.}
  \item{and.split}{The AND-operator (if any).}
  \item{use.tilde}{Logical, use '~' for negation with bivalent variables.}
  \item{x}{An object of class 'DeMorgan'.}
}

\details{
The two De Morgan laws posit that the negation of a disjunction is the 
conjunction of its separate negations, and the negation of a conjunction is the 
disjunction of its separate negations (Hohn 1966, p.80).

The argument \code{expression} can be any complex string representing a 
Boolean expression of disjunctions and conjunctions, or a solution object of 
class 'qca' (objects returned by the '\code{eQMC}' function).
}

\value{A list of solutions with their negations as components if \code{expression} 
       is an object of class 'qca', or simply a list with the following components 
       if \code{expression} is a string:\cr
  \item{initial}{The initial expression.}
  \item{negated}{The negation of the initial expression.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming, testing\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Hohn, Franz E. 1966. \emph{Applied Boolean Algebra: An Elementary Introduction}. 
2nd ed. New York: Macmillan.

Ragin, Charles C. 1987. \emph{The Comparative Method: Moving beyond Qualitative 
and Quantitative Strategies}. Berkeley: University of California Press.
}

\seealso{\code{\link{eQMC}}}

\examples{
# example from Ragin (1987, p.99)
DeMorgan("AC + B~C")

# with different AND-operators
DeMorgan("A*C + B*~C", and.split = "*")
DeMorgan("A&C + B&~C", and.split = "&")

# use solution object of class 'qca' returned by 'eQMC' function, 
# even with multiple models
data(d.represent)
KRO.ps <- eQMC(d.represent, outcome = "WNP")
DeMorgan(KRO.ps)
}

\keyword{functions}
