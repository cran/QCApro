\name{QCApro-package}

\alias{QCApro-package}

\alias{QCApro}

\docType{package}

\title{
Advanced Functionality for Performing and Evaluating Qualitative Comparative Analysis
}

\description{
\pkg{QCApro} is a successor package to the \pkg{QCA} package, with \pkg{QCA} 1.1-4 as its original basis (\href{http://www.alrik-thiem.net/software/}{Dusa and Thiem 2014}; \href{http://www.maxwell.syr.edu/uploadedFiles/moynihan/cqrm/Newsletter10_2.pdf}{Thiem and Dusa 2012}; \href{http://dx.doi.org/10.1177/0894439313478999}{2013a}; \href{https://journal.r-project.org/archive/2013/RJ-2013-009/index.html}{2013b}; \href{http://www.springer.com/gb/book/9781461445838}{2013c}). Just like its predecessor, \pkg{QCApro} implements the method of \dfn{Qualitative Comparative Analysis} (QCA)---a family of techniques for analyzing configurational data in accordance with the INUS theory of causation (Mackie 1965; 1974), but it has fixed various technical and methodological problems of the \pkg{QCA} package and includes many new features and enhancements for applying QCA. 

Moreover, \pkg{QCApro} is currently the only QCA software that provides many purpose-built functions for testing methodological properties of QCA and QCA-related procedures. For example, the effects of changing discretionary parameters such as the inclusion cut-off on the degree of ambiguity affecting a QCA solution can be analyzed (\href{http://dx.doi.org/10.1177/0049124115610351}{Baumgartner and Thiem 2017a}), the consequences of increasing limited empirical diversity on the probability of QCA not committing causal fallacies can be computed (\href{http://dx.doi.org/10.1177/0049124117701487}{Baumgartner and Thiem 2017b}), and the relation between correlational and implicational independence can be examined (\href{http://dx.doi.org/10.1177/0081175016654736}{Thiem and Baumgartner 2016}).  

Three variants can currently be processed by \pkg{QCApro}: \dfn{crisp-set QCA} (csQCA; Ragin 1987), \dfn{multi-value QCA} (mvQCA; Cronqvist and Berg-Schlosser 2009; \href{http://dx.doi.org/10.1177/1525822x13478135}{Thiem 2013}; \href{http://dx.doi.org/10.1007/s11135-014-0015-x}{2014}) and \dfn{fuzzy-set QCA} (fsQCA; Ragin 2000; 2008). A subvariant of csQCA called \dfn{temporal QCA} (tQCA) is also available (\href{http://dx.doi.org/10.1177/0049124105277197}{Caren and Panofsky 2005}; \href{http://dx.doi.org/10.1177/0049124107313903}{Ragin and Strand 2008}).

Several datasets from various areas are integrated in \pkg{QCApro} so as to facilitate familiarization with the package's functionality. Currently covered are business, management and organization (\code{\link{d.stakeholder}}), education (\code{\link{d.education}}), environmental sciences (\code{\link{d.biodiversity}}), evaluation (\code{\link{d.transport}}), legal studies (\code{\link{d.napoleon}}), political science (\code{\link{d.jobsecurity}}, \code{\link{d.partybans}}, \code{\link{d.represent}}), public health (\code{\link{d.health}}, \code{\link{d.tumorscreen}}), urban affairs (\code{\link{d.urban}}), and sociology (\code{\link{d.homeless}}, \code{\link{d.socialsecurity}}). For more details, see the datasets' documentation files.

As an additional resource, \pkg{QCApro} includes a comprehensive glossary for Configurational Comparative Methods. The glossary is directly accessible via the link 'User guides, package vignettes and other documentation' in the package's help index or the 'doc' folder of the package's installation folder.

If you make use of the \pkg{QCApro} package in your work, please acknowledge it in the interest of good scientific practice and transparency. The package citation displays on loading the package or by using the command \code{citation(package = "QCApro")} after loading. The aforesaid command also provides a suitable BibTeX entry. To browse the latest news about the \pkg{QCApro} package (bug fixes, enhancements, etc.), enter \code{news(package = "QCApro")}.

Happy QCAing!
}

\details{
    \tabular{ll}{
    Package: \tab QCApro\cr
    Type:    \tab Package\cr
    Version: \tab 1.1-2\cr
    Date:    \tab 2018-01-10\cr
    License: \tab GPL-3\cr
    }
}

\author{

\bold{Author}:\cr
Alrik Thiem\cr
Department of Political Science\cr
University of Lucerne, Switzerland\cr
\href{mailto:alrik.thiem@unilu.ch}{E-Mail}\cr
\href{http://www.alrik-thiem.net}{Personal Website}\cr
\href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website}

\bold{Maintainer}:\cr
Alrik Thiem\cr
Department of Political Science\cr
University of Lucerne, Switzerland\cr
\href{mailto:alrik.thiem@unilu.ch}{E-Mail}\cr
\href{http://www.alrik-thiem.net}{Personal Website}\cr
\href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website}
}

\references{
Baumgartner, Michael, and Alrik Thiem. 2017a. \dQuote{Model Ambiguities in Configurational Comparative Research.} \emph{Sociological Methods & Research} \strong{46} (4):954-87. DOI: \href{http://dx.doi.org/10.1177/0049124115610351}{10.1177/0049124115610351}.

Baumgartner, Michael, and Alrik Thiem. 2017b. \dQuote{Often Trusted but Never 
(Properly) Tested: Evaluating Qualitative Comparative Analysis.} \emph{Sociological Methods & Research}. Advance online publication. DOI: \href{http://dx.doi.org/10.1177/0049124117701487}{10.1177/0049124117701487}.

Caren, Neal, and Aaron Panofsky. 2005. \dQuote{TQCA: A Technique for Adding 
Temporality to Qualitative Comparative Analysis.} 
\emph{Sociological Methods & Research} \strong{34} (2):147-72.
DOI: \href{http://dx.doi.org/10.1177/0049124105277197}{10.1177/0049124105277197}.

Cronqvist, Lasse, and Dirk Berg-Schlosser. 2009. \dQuote{Multi-Value QCA (mvQCA).} 
In \emph{Configurational Comparative Methods: Qualitative Comparative Analysis (QCA) 
and Related Techniques}, eds. B. Rihoux and C. C. Ragin. London: Sage Publications, 
pp. 69-86.

Dusa, Adrian, and Alrik Thiem. 2014. \emph{QCA: A Package for Qualitative Comparative Analysis}. R Package Version 1.1-4. URL: \url{http://www.alrik-thiem.net/software/}.

Mackie, John L. 1965. \dQuote{Causes and Conditions.} \emph{American Philosophical Quarterly} \strong{2} (4):245-64. URL: \href{http://www.jstor.org/stable/20009173}{http://www.jstor.org/stable/20009173}.

Mackie, John L. 1974. \emph{The Cement of the Universe: A Study of Causation}. Oxford: Oxford University Press.

Ragin, Charles C. 1987. \emph{The Comparative Method: Moving beyond Qualitative and 
Quantitative Strategies}. Berkeley: University of California Press.

Ragin, Charles C. 2000. \emph{Fuzzy-Set Social Science}. Chicago: University of 
Chicago Press.

Ragin, Charles C. 2008. \emph{Redesigning Social Inquiry: Fuzzy Sets and Beyond}. 
Chicago: University of Chicago Press.

Ragin, Charles C., and Sarah Ilene Strand. 2008. \dQuote{Using Qualitative 
Comparative Analysis to Study Causal Order: Comment on Caren and Panofsky (2005).} 
\emph{Sociological Methods & Research} \strong{36} (4):431-41.
DOI: \href{http://dx.doi.org/10.1177/0049124107313903}{10.1177/0049124107313903}.

Thiem, Alrik. 2013. \dQuote{Clearly Crisp, and Not Fuzzy: A Reassessment of the (Putative) Pitfalls of Multi-Value QCA.} \emph{Field Methods} \strong{25} (2):197-207. DOI: \href{http://dx.doi.org/10.1177/1525822x13478135}{10.1177/1525822x13478135}.

Thiem, Alrik. 2014. \dQuote{Parameters of Fit and Intermediate Solutions in 
Multi-Value Qualitative Comparative Analysis.} \emph{Quality & Quantity} \strong{49} (2):657-74. DOI: \href{http://dx.doi.org/10.1007/s11135-014-0015-x}{10.1007/s11135-014-0015-x}.

Thiem, Alrik, and Adrian Dusa. 2012. \dQuote{Introducing the QCA Package: A Market Analysis and Software Review.} \emph{Qualitative & Multi-Method Research} \strong{10} (2):45-9. \href{http://www.maxwell.syr.edu/uploadedFiles/moynihan/cqrm/Newsletter10_2.pdf}{Link}.

Thiem, Alrik, and Adrian Dusa. 2013a. \dQuote{Boolean Minimization in Social 
Science Research: A Review of Current Software for Qualitative Comparative 
Analysis (QCA).} \emph{Social Science Computer Review} \strong{31} (4):505-21.
DOI: \href{http://dx.doi.org/10.1177/0894439313478999}{10.1177/0894439313478999}.

Thiem, Alrik, and Adrian Dusa. 2013b. \dQuote{QCA: A Package for Qualitative 
Comparative Analysis.} \emph{The R Journal} \strong{5} (1):87-97. \href{https://journal.r-project.org/archive/2013/RJ-2013-009/index.html}{Link}.

Thiem, Alrik, and Adrian Dusa. 2013c. \emph{Qualitative Comparative Analysis with R:
A User's Guide}. New York: Springer. \href{http://www.springer.com/gb/book/9781461445838}{Link}.

Thiem, Alrik, and Michael Baumgartner. 2016. \dQuote{Modeling Causal Irrelevance in Evaluations of Configurational Comparative Methods.} \emph{Sociological Methodology} 46 (1):345-57. DOI: \href{http://dx.doi.org/10.1177/0081175016654736}{10.1177/0081175016654736}.
}

\keyword{package}