\name{QCApro-package}

\alias{QCApro-package}

\alias{QCApro}

\docType{package}

\title{
Professional Functionality for Performing and Evaluating Qualitative Comparative Analysis
}

\description{
The \pkg{QCApro} package is the direct successor to the \pkg{QCA} package (\href{http://cran.r-project.org/package=QCA}{Dusa and Thiem 2014}; \href{http://www.maxwell.syr.edu/uploadedFiles/moynihan/cqrm/Newsletter10_2.pdf}{Thiem and Dusa 2012}; \href{http://dx.doi.org/10.1177/0894439313478999}{2013a}; \href{http://journal.r-project.org/archive/2013-1/thiem-dusa.pdf}{2013b}; \href{http://www.springer.com/gb/book/9781461445838}{2013c}). Just like its predecessor, \pkg{QCApro} implements the method of \dfn{Qualitative Comparative Analysis} (QCA)---a family of techniques for analyzing configurational data in accordance with the INUS theory of causation (Mackie 1965; 1974), but it has fixed various 
methodological errors and includes many new features and enhancements for applying and testing QCA.

Three variants can currently be processed by \pkg{QCApro}: \dfn{crisp-set QCA} (csQCA; Ragin 1987), \dfn{multi-value QCA} (mvQCA; Cronqvist and Berg-Schlosser 2009; \href{http://dx.doi.org/10.1177/1525822x13478135}{Thiem 2013}; \href{http://dx.doi.org/10.1007/s11135-014-0015-x}{2014}) and \dfn{fuzzy-set QCA} (fsQCA; Ragin 2000; 2008). A subvariant of csQCA called \dfn{temporal QCA} (tQCA) is also available (\href{http://dx.doi.org/10.1177/0049124105277197}{Caren and Panofsky 2005}; \href{http://dx.doi.org/10.1177/0049124107313903}{Ragin and Strand 2008}).

Several data-sets from various areas are integrated in \pkg{QCApro} so as to facilitate familiarization with the package's functionality. Currently covered are business, management and organization (\code{\link{d.stakeholder}}), education (\code{\link{d.education}}), environmental sciences (\code{\link{d.biodiversity}}), evaluation (\code{\link{d.transport}}), legal studies (\code{\link{d.napoleon}}), political science (\code{\link{d.jobsecurity}}, \code{\link{d.partybans}}, \code{\link{d.represent}}), public health (\code{\link{d.health}}), urban affairs (\code{\link{d.urban}}), and sociology (\code{\link{d.homeless}}, \code{\link{d.socialsecurity}}). For more details, see the data-sets' documentation files. Please also note that many more data-sets for QCA are available on the COMPASSS website at \url{http://www.compasss.org/bibdata.htm}.

If you make use of the \pkg{QCApro} package in your work, please acknowledge it in the interest of good scientific conduct and transparency. The package citation displays on loading the package or by using the command \code{citation(package = "QCApro")} after loading. The aforesaid command also provides a suitable BibTeX entry. To browse the latest news about the \pkg{QCApro} package (bug fixes, enhancements, etc.), enter \code{news(package = "QCApro")}. For version 1.0-0, this will detail all changes in relation to the predecessor package \pkg{QCA} version 1.1-4.

Happy QCAing!
}

\details{
    \tabular{ll}{
    Package: \tab QCApro\cr
    Type:    \tab Package\cr
    Version: \tab 1.0-0\cr
    Date:    \tab 2016-01-27\cr
    License: \tab GPL-3\cr
    }
}

\author{

\bold{Author}:\cr
Alrik Thiem\cr
Department of Philosophy\cr
University of Geneva, Switzerland\cr
\href{mailto:alrik.thiem@unige.ch}{E-Mail}\cr
\href{http://www.alrik-thiem.net}{Personal Website}\cr
\href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website}

\bold{Maintainer}:\cr
Alrik Thiem\cr
Department of Philosophy\cr
University of Geneva, Switzerland\cr
\href{mailto:alrik.thiem@unige.ch}{E-Mail}\cr
\href{http://www.alrik-thiem.net}{Personal Website}\cr
\href{https://www.researchgate.net/profile/Alrik_Thiem}{ResearchGate Website}
}

\references{
Caren, Neal, and Aaron Panofsky. 2005. \dQuote{TQCA: A Technique for Adding 
Temporality to Qualitative Comparative Analysis.} 
\emph{Sociological Methods & Research} \strong{34} (2):147-72.
DOI: \href{http://dx.doi.org/10.1177/0049124105277197}{10.1177/0049124105277197}.

Cronqvist, Lasse, and Dirk Berg-Schlosser. 2009. \dQuote{Multi-Value QCA (mvQCA).} 
In \emph{Configurational Comparative Methods: Qualitative Comparative Analysis (QCA) 
and Related Techniques}, eds. B. Rihoux and C. C. Ragin. London: Sage Publications, 
pp. 69-86.

Dusa, Adrian, and Alrik Thiem. 2014. \emph{QCA: A Package for Qualitative Comparative Analysis}. R Package Version 1.1-4. URL: \href{http://cran.r-project.org/package=QCA}{http://cran.r-project.org/package=QCA}.

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
Comparative Analysis.} \emph{The R Journal} \strong{5} (1):87-97. \href{http://journal.r-project.org/archive/2013-1/thiem-dusa.pdf}{Link}.

Thiem, Alrik, and Adrian Dusa. 2013c. \emph{Qualitative Comparative Analysis with R:
A User's Guide}. New York: Springer. \href{http://www.springer.com/gb/book/9781461445838}{Link}.
}

\keyword{package}