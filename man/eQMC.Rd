\name{eQMC}

\alias{eQMC}
\alias{is.qca}

\title{Minimization with Enhanced Quine-McCluskey Algorithm}

\description{
This function performs the minimization. Although it is called '\emph{e}QMC', the 
implemented algorithm is different from the classical Quine-McCluskey (QMC) 
algorithm. Instead of QMC's approach of using positive min-terms and remainders
to perform minimization, \emph{e}QMC uses positive and negative min-terms, but no 
remainders. See \href{http://dx.doi.org/10.1080/0022250X.2014.897949}{Dusa and Thiem (2015)} and \href{http://dx.doi.org/10.1177/0049124115589032}{Thiem (2015)} for more details.
}

\usage{
eQMC(data, outcome = c(""), neg.out = FALSE, exo.facs = c(""), 
     relation = "suf", n.cut = 1, incl.cut1 = 1, incl.cut0 = 1, 
     minimize = c("1"), sol.type = "ps", row.dom = FALSE, 
     min.dis = FALSE, omit = c(), dir.exp = c(), details = FALSE,
     show.cases = FALSE, inf.test = c(""), use.tilde = FALSE, 
     use.letters = FALSE, ...)

is.qca(x)
}

\arguments{
  \item{data}{A truth table object or a set of configurational data (of class 
       'matrix' or 'data.frame').}
  \item{outcome}{A vector of outcomes.}
  \item{neg.out}{Logical, use negation of \code{outcome} (ignored if \code{data}
        is a truth table object).}
  \item{exo.facs}{A character vector with the names of the exogenous factors.}
  \item{relation}{The required relation of a model antecendent to the 
        \code{outcome}; either \code{"suf"} (only sufficiency required) 
        or \code{"sufnec"} (both sufficiency and necessity required).}
  \item{n.cut}{The minimum number of cases with set membership score above 0.5 
        for an output function value of "0", "1" or "C"; an integer between 1 and
        the maximum number of cases for all non-remainder min-terms.}
  \item{incl.cut1}{The minimum sufficiency inclusion score for an output function 
        value of "1".}
  \item{incl.cut0}{The maximum sufficiency inclusion score for an output function 
        value of "0".}
  \item{minimize}{A vector of output function values for which a solution is sought.}
  \item{sol.type}{A character scalar specifying the QCA solution type that should 
        be applied; either "ps" (parsimonious solution), "ps+" (parsimonious 
        solution including both positive and contradiction min-terms), "cs" (
        conservative solution) or "cs+" (conservative solution including both 
        positive and contradiction min-terms). Note that only "ps" and "ps+" 
        generate correct solutions.}
  \item{row.dom}{Logical, impose row dominance as a constraint on the solution to
        eliminate dominated inessential prime implicants. For causal data analysis,
        this argument must be set to \code{FALSE}.}
  \item{min.dis}{Logical, impose minimal disjunctivity as a constraint on the 
        solution to eliminate models with more prime implicants than the model(s) 
        with the fewest prime implicants. For causal data analysis, this argument 
        must be set to \code{FALSE}.}
  \item{omit}{A vector of min-term index values or a matrix of min-terms to be 
        omitted from minimization.}
  \item{dir.exp}{A vector of directional expectations for deriving intermediate 
        solutions; can only be used in conjunction with \code{sol.type = "ps"}
        or \code{sol.type = "ps+"}. Note that neither conservative nor intermediate
        solutions produce correct solutions. This argument is only retained for
        purposes of method evaluation.}
  \item{details}{Logical, present solution details (inclusion, raw coverage and
        unique coverage scores).}
  \item{show.cases}{Logical, also print case names as part of a solution's details;
        \code{details} must be set to \code{TRUE} (do not use this option with many
        cases and/or long case names).}
  \item{inf.test}{A vector of length two specifying the inference-statistical 
        test to be performed (currently only \code{"binom"}) and the critical
        significance level.}
  \item{use.tilde}{Logical, use tilde operator ("~") for negation with bivalent
        (crisp-set and fuzzy-set) factors.}
  \item{use.letters}{Logical, use single letters (in alphabetical order) instead 
        of original variable names.}
  \item{...}{Other arguments.}
  \item{x}{An object of class 'qca'.}
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
factor (e.g., \code{O} as a short-cut for \code{O{1}}). Outcomes from multivalent crisp-set factors always require curly-bracket notation. Outcomes can be single 
levels of factors not simultaneously passed to \code{exo.facs}, or levels from 
any subset of the factors specified in \code{exo.facs} if \code{data} is not a 
truth table object. At least one outcome always has to be specified.

If multiple outcomes are specified, their factors must also be specified in 
\code{exo.facs}. In this case, solution details will not be printed by default 
(see the example on mimicking Coincidence Analysis below).

The logical argument \code{neg.out} controls whether \code{outcome} is to be 
analyzed or its negation. If \code{outcome} is a level from a multivalent factor, 
\code{neg.out = TRUE} makes the disjunction of all remaining levels the outcome.

The argument \code{exo.facs} specifies the exogenous factors. If omitted, all 
factors in \code{data} are used except that of the \code{outcome}. With multiple outcomes, all factors in \code{data} are used. Please note that computation times 
may increase significantly beyond 17 exogenous factors, and that the computation 
of a solution may not be possible at all depending on end-user machine constraints.

The argument \code{relation} specifies the relation between the antecedent of a
model and the outcome. It accepts either the value \code{"suf"} or \code{"sufnec"}. 
If \code{relation = "suf"} (default), only sufficiency is used as a criterion in 
identifying a model. If \code{relation = "sufnec"}, models must be sufficient and 
necessary for the outcome to be identified. The argument \code{incl.cut1} then 
acts as the cut-off for the sufficiency inclusion of a min-term as well as the 
necessity inclusion of the final model(s).

Min-terms that contain fewer than \code{n.cut} cases with membership scores above 
0.5 are coded as remainders (\code{OUT = "?"}). If the number of such cases is at 
least \code{n.cut}, min-terms with an inclusion score of at least \code{incl.cut1} 
are coded positive (\code{OUT = "1"}), min-terms with an inclusion score below 
\code{incl.cut1} but with at least \code{incl.cut0} are coded as a contradiction 
(\code{OUT = "C"}), and min-terms with an inclusion score below \code{incl.cut0} 
are coded negative (\code{OUT = "0"}). If \code{incl.cut0} is not explicitly 
changed, it is set equal to \code{incl.cut1}.

The argument \code{minimize} specifies a vector of suitable values of the output
function for which a solution is sought. Vectors of such values are \code{"1"}
(default; positive min-terms), \code{"C"} (contradictions), \code{"0"} (negative
min-terms), \code{c("1", "C")} and \code{c("0", "C")}, but not \code{c("1", "0")} 
and \code{c("1", "0", "C")}. Note that for \code{"0"}, \code{"C"} and 
\code{c("0", "C")}, the respective min-terms will be processed but no solution 
details will be printed. Also note that \code{minimize = "0"} is not the same
as using \code{neg.out = TRUE}.

The argument \code{sol.type} specifies the QCA solution type that should be 
applied. It accepts either \code{"ps"} (default, parsimonious solution), 
\code{"ps+"} (parsimonious solution including both positive min-terms and 
contradictions), \code{"cs"} (conservative solution) or \code{"cs+"} (conservative 
solution including both positive min-terms and contradictions). As only the 
parsimonious solution generates methodologically correct solutions ( \href{https://www.researchgate.net/publication/280555425_Often_Trusted_But_Never_Properly_Tested_Evaluating_Qualitative_Comparative_Analysis}{Baumgartner 
and Thiem 2015a}), \code{sol.type} should normally be left unchanged.

The logical argument \code{row.dom} controls whether the principle of row dominance
is imposed as a constraint on the solution. An inessential prime implicant \eqn{P}
dominates another \eqn{Q} if all configurations covered by \eqn{Q} are also
covered by \eqn{P}, but they are not interchangeable (cf. \href{http://dx.doi.org/10.1002/j.1538-7305.1956.tb03835.x}{McCluskey 1956, 1425};
McCluskey 1965, 164-152). If row dominance is operative, models that contain
dominated prime implicants will not be returned. For purposes of causal data
analysis, \code{row.dom} must be set to \code{FALSE}.

The logical argument \code{min.dis} controls whether the principle of minimal
disjunctivity is imposed as a constraint on the solution (McCluskey 1965, 123-126).
If minimal disjunctivity is operative, models that contain more than the number
of prime implicants of the model(s) with the fewest prime implicants will not be
returned. \emph{For purposes of causal data analysis, both} \code{row.dom} 
\emph{and} \code{min.dis} \emph{must be set to} \code{FALSE} (\href{http://dx.doi.org/10.1177/0049124115610351}{Baumgartner and Thiem 2015b}; \href{http://dx.doi.org/10.1177/0193841x14550863}{Thiem 2014b}).

The argument \code{omit} can be used to omit min-terms from the minimization 
process \emph{ex ante}. It accepts a vector of row numbers from the truth table 
or a matrix of min-terms of the same order as passed to the \code{truthTable} 
function (if the argument \code{data} is a truth table object) or as specified
in the argument \code{exo.facs}.

Neither conservative nor intermediate solutions produce correct solutions (\href{https://www.researchgate.net/publication/280555425_Often_Trusted_But_Never_Properly_Tested_Evaluating_Qualitative_Comparative_Analysis}{Baumgartner and Thiem 2015a}). The 
\code{dir.exp} argument is retained only for purposes of method evaluation in 
relation to intermediate solutions. It specifies directional expectations for 
separating easy from difficult counterfactuals in simplifying assumptions. For 
bivalent crisp and fuzzy-set factors, expectations should be specified as a vector 
of the same length and the same order of condition variables as provided in 
\code{exo.facs}. For bivalent factors, a value of either "0" or "1" indicates
that the corresponding factor is expected to contribute to a positive output
function value, while a dash, "-", indicates that one or the other level of the
corresponding factor does so. For multivalent factors, multiple levels have to 
be enclosed by double quotes and separated by a semicolon (see mvQCA example using 
\href{http://dx.doi.org/10.1080/13510347.2010.491189}{Hartmann and Kemmerzell (2010)} 
below). In some situations, directional expectations in mvQCA generate easy 
counterfactuals that do not contribute to parsimony 
(\href{http://dx.doi.org/10.1007/s11135-014-0015-x}{Thiem 2014a}).

If \code{details = TRUE}, parameters of fit (inclusion, raw coverage, and unique 
coverage) will be printed for each solution and its respective prime implicants. 
Essential prime implicants are listed first in the solution output and in the top 
part of the parameters-of-fit table. Inessential prime implicants are listed in 
brackets in the solution output and in the middle part of the parameters-of-fit 
table, together with their unique coverage scores under each individual model.
Inclusion and coverage scores for each model are provided in the bottom part
of the parameters-of-fit table.

The logical argument \code{show.cases} controls whether case names are displayed
next to their corresponding prime implicants (do not use with many cases and/or 
long case names!). In the parameters-of-fit table, semicolons separate cases from 
different min-terms, whereas commas separate cases from the same min-term.

The argument \code{inf.test} provides functionality for basing output function
value codings on inference-statistical tests. Currently, only an exact binomial
test (\code{"binom"}) is available, which requires the data to contain only
bivalent or multivalent crisp-set factors. The argument requires a vector of
length two, comprising the test and a critical significance level. If the
empirical inclusion score of a min-term is not significantly lower than 
\code{incl.cut1}, it will be coded positive (\code{OUT = "1"}). If it is 
significantly lower than \code{incl.cut1} yet still significantly higher than 
\code{incl.cut0}, it will be coded as a contradiction (\code{OUT = "C"}). If it
is not significantly higher than \code{incl.cut0}, it will be coded negative
(\code{OUT = "0"}).

The argument \code{use.tilde} should only be used for bivalent factors. If the exogenous factors are already named with single letters, the argument 
\code{use.letters} will have no effect when set to \code{TRUE}. Otherwise, 
upper-case letters will replace original factor names in alphabetical order.
}

\value{An object of class 'qca' for single outcomes and 'mqca' for multiple
       outcomes. Objects of class 'qca' are lists with the following 
       ten main components:\cr
  \item{tt}{The truth table object.}
  \item{excluded}{The line numbers of the negative min-terms.}
  \item{initials}{The positive (non-remainder) min-terms.}
  \item{PIs}{The prime implicants.}
  \item{PIchart}{The list of prime implicant charts.}
  \item{solution}{The list of solutions.}
  \item{essential}{The list of essential prime implicants.}
  \item{pims}{The list of model prime implicant set membership scores.}
  \item{SA}{The list of simplifying assumptions that would have been used by 
        Quine-McCluskey minimization.}
  \item{i.sol}{A list of components specific to intermediate solution(s), 
        including the prime implicant chart, model prime implicant membership 
        scores, (non-simplifying) easy counterfactuals and difficult 
        counterfactuals.}
}

\author{
Alrik Thiem (\href{http://www.alrik-thiem.net}{Personal Website}; \href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website})
}

\section{Contributors}{\tabular{ll}{
Dusa, Adrian\tab: development, programming\cr
Thiem, Alrik\tab: development, documentation, testing
}}

\references{
Baumgartner, Michael. 2009. \dQuote{Inferring Causal Complexity.} 
\emph{Sociological Methods & Research} \strong{38} (1):71-101.
DOI: \href{http://dx.doi.org/10.1177/0049124109339369}{10.1177/0049124109339369}.

Baumgartner, Michael, and Alrik Thiem. 2015a. \emph{Often Trusted but Never 
(Properly) Tested: Evaluating Qualitative Comparative Analysis}. Paper presented 
at the 12th Conference of the European Sociological Association, 25-28 August, 
Czech Technical University, Prague (Czech Republic). \href{https://www.researchgate.net/publication/280555425_Often_Trusted_But_Never_Properly_Tested_Evaluating_Qualitative_Comparative_Analysis}{Link}.

Baumgartner, Michael, and Alrik Thiem. 2015b. \dQuote{Model Ambiguities in Configurational Comparative Research.} \emph{Sociological Methods & Research}. Advance online publication. DOI: \href{http://dx.doi.org/10.1177/0049124115610351}{10.1177/0049124115610351}.

Dusa, Adrian. 2007. \dQuote{User Manual for the QCA(GUI) Package in R.} \emph{Journal of Business Research} \strong{60} (5):576-86. DOI: \href{http://dx.doi.org/10.1016/j.jbusres.2007.01.002}{10.1016/j.jbusres.2007.01.002}.

Dusa, Adrian, and Alrik Thiem. 2015. \dQuote{Enhancing the Minimization of Boolean 
and Multivalue Output Functions with eQMC.} \emph{Journal of Mathematical Sociology} 
\strong{39} (2):92-108. DOI: 
\href{http://dx.doi.org/10.1080/0022250X.2014.897949}{10.1080/0022250X.2014.897949}.

Emmenegger, Patrick. 2011. \dQuote{Job Security Regulations in Western Democracies: 
A Fuzzy Set Analysis.} \emph{European Journal of Political Research} \strong{50} (3):336-64.
DOI: \href{http://dx.doi.org/10.1111/j.1475-6765.2010.01933.x}{10.1111/j.1475-6765.2010.01933.x}.

Hartmann, Christof, and Joerg Kemmerzell. 2010. \dQuote{Understanding Variations 
in Party Bans in Africa.} \emph{Democratization} \strong{17} (4):642-65.
DOI: \href{http://dx.doi.org/10.1080/13510347.2010.491189}{10.1080/13510347.2010.491189}.

Huang, Ronggui. 2014. \emph{QCA3: Yet another Package for Qualitative Comparative Analysis}. R Package Version 0.0-7. URL: \href{http://cran.r-project.org/package=QCA3}{http://cran.r-project.org/package=QCA3}.

McCluskey, Edward J. 1956. \dQuote{Minimization of Boolean Functions.} 
\emph{Bell Systems Technical Journal} \strong{35} (6):1417-44. DOI: \href{http://dx.doi.org/10.1002/j.1538-7305.1956.tb03835.x}{10.1002/j.1538-7305.1956.tb03835.x}.

McCluskey, Edward J. 1965. \emph{Introduction to the Theory of Switching Circuits}.
Princeton: Princeton University Press.

Krook, Mona Lena. 2010. \dQuote{Women's Representation in Parliament: A Qualitative 
Comparative Analysis.} \emph{Political Studies} \strong{58} (5):886-908.
DOI: \href{http://dx.doi.org/10.1111/j.1467-9248.2010.00833.x}{10.1111/j.1467-9248.2010.00833.x}.

Ragin, Charles C. 2008. \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}. Chicago: University of Chicago Press. \href{http://press.uchicago.edu/ucp/books/book/chicago/R/bo5973952.html}{Link}.

Schneider, Carsten Q., and Claudius Wagemann. 2012. \emph{Set-Theoretic Methods for the Social Sciences: A Guide to Qualitative Comparative Analysis (QCA)}. Cambridge: Cambridge University Press. \href{http://www.cambridge.org/ch/academic/subjects/politics-international-relations/research-methods-politics/set-theoretic-methods-social-sciences-guide-qualitative-comparative-analysis?format=PB}{Link}.

Thiem, Alrik. 2014a. \dQuote{Parameters of Fit and Intermediate Solutions in Multi-Value
Qualitative Comparative Analysis.} \emph{Quality & Quantity} \strong{49} (2):657-74.
DOI: \href{http://dx.doi.org/10.1007/s11135-014-0015-x}{10.1007/s11135-014-0015-x}.

Thiem, Alrik. 2014b. \dQuote{Navigating the Complexities of Qualitative Comparative 
Analysis: Case Numbers, Necessity Relations, and Model Ambiguities.} \emph{Evaluation Review} \strong{38} (6):487-513. 
DOI: \href{http://dx.doi.org/10.1177/0193841x14550863}{10.1177/0193841x14550863}.

Thiem, Alrik. 2015. \dQuote{Using Qualitative Comparative Analysis for Identifying 
Causal Chains in Configurational Data: A Methodological Commentary on Baumgartner 
and Epple (2014).} \emph{Sociological Methods & Research} \strong{44} (4):723-36.
DOI: \href{http://dx.doi.org/10.1177/0049124115589032}{10.1177/0049124115589032}.
}

\seealso{\code{\link{pof}}, \code{\link{truthTable}}}

\examples{
# csQCA using Krook (2010)
#-------------------------
data(d.represent)
head(d.represent)

# solution with details and case names
KRO <- eQMC(d.represent, outcome = "WNP", details = TRUE, show.cases = TRUE)
KRO

# check PI chart
KRO$PIchart

# solution with truth table object
KRO.tt <- truthTable(d.represent, outcome = "WNP")
KRO <- eQMC(KRO.tt)
KRO

# simplifying assumptions (SAs) that would have been used with Quine-McCluskey 
# optimization
KRO$SA

# fsQCA using Emmenegger (2011)
#------------------------------
data(d.jobsecurity)
head(d.jobsecurity)

# solution with details
EMM <- eQMC(d.jobsecurity, outcome = "JSR", incl.cut1 = 0.9, details = TRUE)
EMM

# are the model prime implicants also sufficient for the negation of the outcome?
pof(EMM$pims, outcome = "JSR", d.jobsecurity, neg.out = TRUE, relation = "suf")

# are the negations of the model prime implicants also sufficient for the outcome?
pof(1 - EMM$pims, outcome = "JSR", d.jobsecurity, relation = "suf")

# plot all three prime implicants of the solution
PIsc <- EMM$pims
par(mfrow = c(1, 3))
for(i in 1:3){
 plot(PIsc[, i], d.jobsecurity$JSR, pch = 19, ylab = "JSR",
  xlab = names(PIsc)[i], xlim = c(0, 1), ylim = c(0, 1),
  main = paste("Prime Implicant", print(i)))
 mtext(paste(
  "Inclusion = ", round(EMM$IC$overall$incl.cov$incl[i], 3),
  "; Coverage = ", round(EMM$IC$overall$incl.cov$cov.r[i], 3)), 
  cex = 0.7, line = 0.4)
 abline(h = 0.5, lty = 2, col = gray(0.5))
 abline(v = 0.5, lty = 2, col = gray(0.5))
 abline(0, 1)
}

# mvQCA using Hartmann and Kemmerzell (2010)
#-------------------------------------------
data(d.partybans)
head(d.partybans)

# specify exogenous factors beforehand
exo.facs <- c("C", "F", "T", "V")

# parsimonious solution with contradictions included
HK.sol <- eQMC(d.partybans, outcome = "PB{1}", exo.facs = exo.facs,
  incl.cut0 = 0.4, sol.type = "ps+", details = TRUE)
HK.sol

# which are the two countries in T{2} but not PB{1}?
rownames(d.partybans[d.partybans$T == 2 & d.partybans$PB != 1, ])

# QCA with multiple outcomes from multivalent variables
#------------------------------------------------------
d.mmv <- data.frame(A = c(2,0,0,1,1,1,2,2), B = c(2,2,2,2,1,1,0,0), 
                    C = c(0,1,0,0,0,2,1,0), D = c(2,1,2,2,3,1,3,0), 
                    E = c(3,2,3,3,0,1,3,2), 
  row.names = letters[1:8])
head(d.mmv)

mmv.s <- eQMC(d.mmv, outcome = c("D{2}", "E{3}"))
mmv.s

# use quotes with curly-bracket notation to access solution component
print(mmv.s$"E{3}", details = TRUE, show.cases = TRUE)

# negation of outcome from multivalent factor is disjunction of all other
# levels; high under-determination (18 models)
mmv.s <- eQMC(d.mmv, outcome = "E{3}", neg.out = TRUE)
mmv.s

# causal chains with QCA (Thiem 2015); data from Baumgartner (2009);
# note that CNA and QCA results are not always equal because CNA applies a
# different concept of the truth table that does not take each min-term's
# individual inclusion score into consideration
#-----------------------------------------------------------------------------
d.Bau <- data.frame(
  U = c(1,1,1,1,0,0,0,0), D = c(1,1,0,0,1,1,0,0),
  L = c(1,1,1,1,1,1,0,0), G = c(1,0,1,0,1,0,1,0),
  E = c(1,1,1,1,1,1,1,0),
  row.names = letters[1:8])
head(d.Bau)

# with multiple outcomes, no solution details are printed;
# "causal-chain structure": (D + U <=> L) * (G + L <=> E)
# "common-cause structure": (D + U <=> L) * (G + D + U <=> E)
Bau.cna <- eQMC(d.Bau, outcome = names(d.Bau), relation = "sufnec")
Bau.cna

# get the truth table, solution details and case names for outcome "E"
print(Bau.cna$E, details = TRUE, show.cases = TRUE)

# examples relating to QCA method evaluation
#-------------------------------------------
#
# is the conservative solution (QCA-CS) really "conservative"?
#-------------------------------------------------------------
# Ragin (2008, 173): "The complex [conservative] solution [...] does not 
# permit any counterfactual cases and thus no simplifying assumptions 
# regarding combinations of conditions that do not exist in the data.";
# the conservative solution is "[c]onservative because [...] the
# researcher [...] is exclusively guided by the empirical information 
# at hand" (Schneider and Wagemann 2012, 162)
#
# in fact, QCA-CS makes extremely strong assumptions on ALL remainders;
# QCA-CS assumes every remainder exists at least 'freq.cut' times, 
# and occurs with the negation of the outcome more than 
# 'freq.cut' * (1 - 'incl.cut1') times

# create a test data-set 'CS' with 32 cases and randomly assign values 
# on the endogenous factor 'Z'
CS <- data.frame(mintermMatrix(rep(2,5)))
CS$Z <- sample(0:1, 2^5, replace = TRUE)

# randomly draw 20 cases to create a limitedly diverse data-set 'CS.LD'
# and turn all 12 remainder min-terms into observations that occur with 
# 'Z = 0' in original data-set 'CS'
CS.LD <- CS[sample(1:2^5, 20), ]
change <- as.numeric(setdiff(rownames(CS), rownames(CS.LD)))
CS$Z[change] <- 0

# create the (conservative) solutions for 'CS' and 'CS.LD'
CS.sol <- eQMC(CS, outcome = "Z")
CS.LD.sol <- eQMC(CS.LD, outcome = "Z", sol.type = "cs")

# test whether the two solutions are identical
identical(unlist(CS.LD.sol$solution), unlist(CS.sol$solution))

# both solutions are identical, for two very different data sets, 
# indicating that QCA-CS draws causal inferences beyond what the data 
# warrants; the lower the diversity index (ratio of non-remainder
# min-terms to all min-terms), the stronger the assumptions QCA-CS makes

\dontrun{
# R packages QCAGUI 2.0 (Dusa 2007) and QCA3 (Huang 2014) 
# produce incorrect solutions
#--------------------------------------------------------
dat <- as.data.frame(mintermMatrix(rep(2, 5)))
dimnames(dat) <- list(as.character(1:2^5), LETTERS[1:5])

# assume data-generating structure is A + B + D <=> E, and mask 
# it by ambiguity
dat <- dat[((dat$A == 1 | dat$B == 1 | dat$D == 1) == dat$E) & 
           ((dat$C == 1 | dat$D == 1) == dat$E), ]

# QCAGUI/QCA 2.0 returns incorrect solutions AC + BC + D <=> E 
# (default settings; onservative solution) and C + D <=> E 
# (non-default settings; parsimonious solution)

install.packages("QCAGUI") # if not yet installed
library(package = "QCAGUI")
QCAGUI::eqmcc(dat, outcome = "E")
QCAGUI::eqmcc(dat, outcome = "E", include = "?")
remove.packages("QCAGUI") # to deinstall QCAGUI

# close R to cleanly unload QCAGUI package components 
# (otherwise conflicts with QCApro may occur)

# QCA3 returns the same incorrect solutions
install.packages("QCA3") # if not yet installed
library(package = "QCA3")

QCA3::reduce(cs_truthTable(dat, outcome = "E", 
             condition = c("A", "B", "C", "D")))
QCA3::reduce(cs_truthTable(dat, outcome = "E", 
             condition = c("A", "B", "C", "D")), remainders = "include")

remove.packages("QCA3") # to deinstall QCA3



}
}

\keyword{functions}
