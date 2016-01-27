\name{NEWS}
\alias{NEWS}

\title{Package News}

\section{Version 1.0-0}{

\subsection{RENAMED}{
\itemize{
\item{function 'allExpressions' to 'implicantMatrix', 'createMatrix' to 
      'mintermMatrix', 'deMorgan' to 'DeMorgan', 'eqmcc' to 'eQMC'}
\item{argument 'conditions' in functions 'eQMC', 'truthTable', and 
      'superSubset' renamed to 'exo.facs' (exogenous factors)}
\item{argument 'include' in function 'eQMC' to 'sol.type' (solution type)}
\item{argument 'explain' in function 'eQMC' to 'minimize'}
\item{argument 'prod.split' in functions 'DeMorgan' and 'factorize' to
      'and.split'}
}}

\subsection{CHANGED}{
\itemize{
\item{functions 'demoChart', 'findSubsets', 'findSupersets', 'getRow', 
      and 'solveChart' removed from user level (internalized)}
\item{error messages corrected and extended in verification functions 
      'verify.data', 'verify.tt', 'verify.qca', 'verify.dir.exp', 
      'verify.mqca', and 'verify.inf.test'}
\item{all documentation files revised and updated}
\item{default value for argument 'min.dis' in function 'eQMC' set to 'TRUE'}
\item{default value for argument 'sol.type' in function 'eQMC' set to 
      'ps' (parsimonious solution)}
\item{truth table not printed any more when dataset is passed directly 
      to function 'eQMC' and argument 'details' is set to 'TRUE'}
}}

\subsection{FIXED}{
\itemize{
\item{bug in functions 'truthTable' and 'eqmcc' (new 'eQMC') produced
      incorrect truth table for outcomes from multivalent factors when
      these outcomes included levels with two or more digits}
\item{bug in function 'eqmcc' ('eQMC') in relation to intermediate solutions}
\item{argument 'inf.test' in function 'truthTable' caused error 
      when argument 'complete = TRUE'}
\item{argument 'data' in function 'eqmcc' (new 'eQMC') caused error when 
      'data' was of class 'matrix' instead of 'data.frame' or 'tt'}
\item{error messages for argument 'data' with respect to missing values 
      restricted to data under factor frame (argument 'exo.facs' plus 
      factor of argument 'outcome') instead of all data under argument 
      'data' (caused error message)}
\item{function 'eqmcc' (new 'eQMC') used original factor labels when
      passed an object of class 'tt' although argument 'use.letters'
      was set to 'TRUE' in call to function 'truthTable'}
\item{function 'print.qca' used incorrect ordering of models in multi-
      model solutions when argument 'relation' was set to 'sufnec' in call
      to function 'eqmcc' (new 'eQMC')}
}}

\subsection{NEW}{
\itemize{
\item{functions 'limitedDiversity', 'implicIndep', 'retention', 'submodels' 
      and 'testTESA' for purposes of method evaluation}
\item{function 'submodels.loop' (internal)}
\item{information in documentation files on contributors (categories: 
      development, documentation, programming, testing)}
}}
}

