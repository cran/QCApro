\name{implicIndep}
\alias{implicIndep}

\title{Test for Implicational Independence between two Factors}

\description{
This evaluation function tests for the implicational independence between two factors. It has been programmed for \href{https://www.researchgate.net/publication/291345976_Modeling_Causal_Irrelevance_in_Evaluations_of_Configurational_Comparative_Methods}{Thiem and Baumgartner (2015)}.
}

\usage{
implicIndep(expression, n.samples = 1, size.sample = 100, corr = "0")
}

\arguments{
  \item{expression}{A string representing the Boolean function to be evaluated.}
  \item{n.samples}{The number of datasets to be sampled.}
  \item{size.sample}{The size of each data sample.}
  \item{corr}{The direction of correlation between the endogenous factor and the
              implicationally independent factor.}
}

\details{
Randomly sample \code{n.samples} different data-sets with uniform probability mass function (any other discrete function would do as well; proficient users may adjust this at the relevant places in the function); and run QCA for each dataset; check whether the irrelevant factor is eliminated and get the correlations between the irrelevant factor and the outcome factor.

The correlation can be controlled with \code{corr}: \code{"0"} means no correlation; \code{"+"} positive correlation, and \code{"-"} negative correlation. The larger the sample size, the larger the positive / negative correlation. The argument \code{expression} may represent any Boolean function in disjunctive normal form as shown below, including proper causal structures such as \code{"(X1*X2 + X3*X4 <=> Y)*(Z1 + z1)"} or non-causal structures such as \code{"(x1*X2 + X1*x2 + X1*X3 + X2*X3 <=> Y)*(Z1 + z1)"} (contains redundant prime implicants) or \code{"(X3*x2 + X2 <=> Y)*(Z1 + z1)"} (contains redundant conjuncts).

If \code{expression} is no causal structure, an additional note will be issued together with the test output for whether the irrelevant factor has been eliminated.

You can use the following possibilities for \code{expression}: \code{"(X1X2 + X3X4 <=> Y)(Z1 + z1)"} if a factor has one letter and a one-digit number, \code{"(AB + CD <=> E)(F + f)"} if a factor has one letter, \code{"(X1*X2 + X3*X4 <=> Y)*(Z1 + z1)"} or \code{"(A*B + C*D <=> E)*(F + f)"} (curly-bracket notation is not supported). Empty spaces and the type of the biconditional operator (\code{<->}/\code{<=>}) have no effect.
}

\value{A list with the following five components:\cr
  \item{tt}{The evaluated truth tables.}
  \item{dat.list}{The test datasets generated on the basis of the evaluaetd truth tables.}
  \item{sol.list}{The corresponding solutions.}
  \item{cor.list}{The correlations between the endogenous factor and the irrelevant factor.}
  \item{test}{The result for whether the irrelevant factor has been eliminated in all
              tests.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Baumgartner, Michael\tab: testing\cr
Thiem, Alrik        \tab: development, documentation, programming, testing
}}

\references{
Thiem, Alrik, and Michael Baumgartner. 2015. \emph{Modelling Causal Irrelevance in Evaluations of Configurational Comparative Methods}. Geneva: Manuscript. \href{https://www.researchgate.net/publication/291345976_Modeling_Causal_Irrelevance_in_Evaluations_of_Configurational_Comparative_Methods}{Link}.
}

\seealso{\code{\link{eQMC}}}

\examples{
\dontrun{
# simulation with 10 sample datasets
simulation <- vector(mode = "list", length = 3)
n.samples <- 10 

# directions of correlation and number of cases in each sample
corr <- c("-", "0", "+") 
nofc <- c(40, 160, 640) 
simulation <- lapply(nofc, function (x) {lapply(corr, function (y) {
  implicIndep("(X1*X2 + X3*X4 <=> Y)*(Z1 + z1)", n.samples, x, y)})})

# has Z1 been eliminated in all data experiments of a block of tests?
series.test <- matrix(sapply(1:length(corr), function (x) {sapply(1:length(nofc),
  function (y) {simulation[[x]][[y]]$test == "Z1 has been eliminated."})}),
  ncol = length(corr), dimnames = list(as.character(nofc), corr))
series.test
}
}

\keyword{functions}


