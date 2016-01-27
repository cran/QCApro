\name{calibrate}

\alias{calibrate}

\title{Calibrate Raw Data into Configurational Data}

\description{
This function generates configurational data from raw data (base variables) and 
some specified threshold(s). The calibration of bivalent fuzzy-set factors is 
possible for positive and negative end-point and mid-point concepts, using the 
method of transformational assignment.
}

\usage{
calibrate(x, type = "crisp", thresholds = NA, include = TRUE, 
          logistic = FALSE, idm = 0.95, ecdf = FALSE, p = 1, q = 1)
}

\arguments{
  \item{x}{An interval or ratio-scaled base variable.}
  \item{type}{The calibration type, either "crisp" or "fuzzy".}
  \item{thresholds}{A vector of thresholds.}
  \item{include}{Logical, include threshold(s) (\code{type = "crisp"} only).}
  \item{logistic}{Calibrate to fuzzy-set variable using the logistic function.}
  \item{idm}{The set inclusion degree of membership for the logistic function.}
  \item{ecdf}{Calibrate to fuzzy-set variable using the empirical cumulative 
        distribution function of the base variable.}
  \item{p}{Parameter: if \eqn{p > 1} concentration, if \eqn{0 < p < 1} dilation 
           below crossover.}
  \item{q}{Parameter: if \eqn{q > 1} dilation, if \eqn{0 < q < 1} concentration 
           above crossover.}
}

\details{
Calibration is the process by which configurational data is produced, that is, 
by which set membership scores are assigned to cases. With interval and ratio-scaled 
base variables, calibration can be based on transformational assignments using 
(piecewise-defined) membership functions.

For \code{type = "crisp"}, one threshold produces a factor with two levels: 0 and 
1. More thresholds produce factors with multiple levels. For example, two thresholds 
produce three levels: 0, 1 and 2.

For \code{type = "fuzzy"}, this function can generate bivalent fuzzy-set variables 
by linear, s-shaped, inverted s-shaped and logistic transformation for end-point 
concepts. It can generate bivalent fuzzy-set variables by trapezoidal, triangular 
and bell-shaped transformation for mid-point concepts (\href{http://www.worldscientific.com/worldscibooks/10.1142/6451}{Bojadziev and Bojadziev 
2007}; \href{http://www.springer.com/de/book/9783540774600}{Clark \emph{et al}. 2008}; \href{http://dx.doi.org/10.1080/13645579.2013.806118}{Thiem 2014}; \href{http://www.springer.com/gb/book/9781461445838}{Thiem and Dusa 2013}).

For calibrating bivalent fuzzy-set variables based on end-point concepts, 
\code{thresholds} should be specified as a numeric vector \code{c(thEX, thCR, thIN)}, 
where \code{thEX} is the threshold for full exclusion, \code{thCR} the threshold 
for the crossover, and \code{thIN} the threshold for full inclusion.

If \code{thEX} \eqn{<} \code{thCR} \eqn{<} \code{thIN}, then the membership 
function is increasing from \code{thEX} to \code{thIN}. If \code{thIN} \eqn{<} 
\code{thCR} \eqn{<} \code{thEX}, then the membership function is decreasing from 
\code{thIN} to \code{thEX}.

For calibrating bivalent fuzzy-set variables based on mid-point concepts, 
\code{thresholds} should be specified as a numeric vector 
\code{c(thEX1, thCR1, thIN1, thIN2, thCR2, thEX2)}, where \code{thEX1} is the first 
(left) threshold for full exclusion, \code{thCR1} the first (left) threshold for 
the crossover, \code{thIN1} the first (left) threshold for full inclusion, 
\code{thIN2} the second (right) threshold for full inclusion, \code{thCR2} the 
second (right) threshold for crossover, and \code{thEX2} the second (right) threshold 
for full exclusion.

If \code{thEX1} \eqn{<} \code{thCR1} \eqn{<} \code{thIN1} \eqn{\le} \code{thIN2}
\eqn{<} \code{thCR2} \eqn{<} \code{thEX2}, then the membership function is first
increasing from \code{thEX1} to \code{thIN1}, then flat between \code{thIN1} and 
\code{thIN2}, and finally decreasing from \code{thIN2} to \code{thEX2}. In contrast, 
if \code{thIN1} \eqn{<} \code{thCR1} \eqn{<} \code{thEX1} \eqn{\le} \code{thEX2}
\eqn{<} \code{thCR2} \eqn{<} \code{thIN2}, then the membership function is first 
decreasing from \code{thIN1} to \code{thEX1}, then flat between \code{thEX1} and 
\code{thEX2}, and finally increasing from \code{thEX2} to \code{thIN2}.

The parameters \code{p} and \code{q} control the degree of concentration and 
dilation. They should be left at their default values unless good reasons for 
changing them exist.

If \code{logistic = TRUE}, the argument \code{idm} specifies the inclusion degree 
of membership.

If \code{ecdf = TRUE}, calibration is based on the empirical cumulative distribution 
function of \code{x}.
}

\value{
A numeric vector of set membership scores between 0 and 1 for bivalent crisp-set 
factors and bivalent fuzzy-set variables, or a numeric vector of levels for 
multivalent crisp-set factors (beginning with 0 at increments of 1).
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: programming\cr
Thiem, Alrik\tab: development, documentation, programming, testing
}}

\references{
Bojadziev, George, and Maria Bojadziev. 2007. \emph{Fuzzy Logic for Business, 
Finance, and Management}. 2nd ed. Hackensack, NJ: World Scientific. \href{http://www.worldscientific.com/worldscibooks/10.1142/6451}{Link}.

Clark, Terry D., Jennifer M. Larson, John N. Mordeson, Joshua D. Potter, and 
Mark J. Wierman. 2008. \emph{Applying Fuzzy Mathematics to Formal Models in 
Comparative Politics}. Berlin: Springer. \href{http://www.springer.com/de/book/9783540774600}{Link}.

Thiem, Alrik. 2014. \dQuote{Membership Function Sensitivity of Descriptive 
Statistics in Fuzzy-Set Relations.} \emph{International Journal of Social 
Research Methodology} \strong{17} (6):625-42. 
DOI: \href{http://dx.doi.org/10.1080/13645579.2013.806118}{10.1080/13645579.2013.806118}.

Thiem, Alrik, and Adrian Dusa. 2013. \emph{Qualitative Comparative Analysis 
with R: A User's Guide}. New York: Springer. \href{http://www.springer.com/gb/book/9781461445838}{Link}.
}

\examples{
# base variable; random draw from standard normal distribution
set.seed(30)
x <- rnorm(30)

# calibration thresholds
th <- quantile(x, seq(from = 0.05, to = 0.95, length = 6))

# calibration of bivalent crisp-set factor
calibrate(x, thresholds = th[3])

# calibration of trivalent crisp-set factor
calibrate(x, thresholds = c(th[2], th[4]))

# fuzzy-set calibration
# 1. positive end-point concept, linear
# 2. positive and corresponding negative end-point concept, logistic
# 3. positive end-point concept, ECDF
# 4. negative end-point concept, s-shaped (quadratic)
# 5. negative end-point concept, inverted s-shaped (root)
# 6. positive mid-point concept, triangular
# 7. positive mid-point concept, trapezoidal
# 8. negative mid-point concept, bell-shaped

yl <- "Set Membership"
xl <- "Base Variable Value"

par(mfrow = c(2,4), cex.main = 1)

plot(x, calibrate(x, type = "fuzzy", thresholds = c(th[1], (th[3]+th[4])/2, 
   th[6])), xlab = xl, ylab = yl, 
   main = "1. positive end-point concept,\nlinear")

plot(x, calibrate(x, type = "fuzzy", thresholds = c(th[1], (th[3]+th[4])/2, 
  th[6]), logistic = TRUE, idm = 0.99), xlab = xl, ylab = yl, 
  main = "2. positive and corresponding negative\nend-point concept, logistic")
  points(x, calibrate(x, type = "fuzzy", thresholds = c(th[6], (th[3]+th[4])/2, 
    th[1]), logistic = TRUE, idm = 0.99))

plot(x, calibrate(x, type = "fuzzy", thresholds = c(th[1], (th[3]+th[4])/2, 
  th[6]), ecdf = TRUE), xlab = xl, ylab = yl, 
  main = "3. positive end-point concept,\nECDF")

plot(x, calibrate(x, type = "fuzzy", thresholds = c(th[6], (th[3]+th[4])/2, 
  th[1]), p = 2, q = 2), xlab = xl, ylab = yl, 
  main = "4. negative end-point concept,\ns-shaped (quadratic)")

plot(x, calibrate(x, type = "fuzzy", thresholds = c(th[6], (th[3]+th[4])/2, 
  th[1]), p = 0.5, q = 0.5), xlab = xl, ylab = yl, 
  main = "5. negative end-point concept,\ninverted s-shaped (root)")

plot(x, calibrate(x, type = "fuzzy", thresholds = th[c(1,2,3,3,4,5)]),
  xlab = xl, ylab = yl, main = "6. positive mid-point concept,\ntriangular")

plot(x, calibrate(x, type = "fuzzy", thresholds = th[c(1,2,3,4,5,6)]),
  xlab = xl, ylab = yl, main = "7. positive mid-point concept,\ntrapezoidal")

plot(x, calibrate(x, type = "fuzzy", thresholds = th[c(3,2,1,5,4,3)],
  p = 3, q = 3), xlab = xl, ylab = yl, 
  main = "8. negative mid-point concept,\nbell-shaped")
}

\keyword{functions}
