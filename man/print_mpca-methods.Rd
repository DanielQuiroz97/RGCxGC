% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\docType{methods}
\name{print_mpca}
\alias{print_mpca}
\alias{print_mpca,MPCA-method}
\title{Print MPCA summary}
\usage{
print_mpca(Object)

\S4method{print_mpca}{MPCA}(Object)
}
\arguments{
\item{Object}{a MPCA object}
}
\description{
`print_mpca` call the MPCA object to print the summary of this
 analysis.
}
\details{
The plot function employs the built-in print function and a
 precomputed MPCA summary to display the explained and cumulative variance
 for each principal component.
}
\examples{

data(MTBLS579)
MTBLS_mpca <- m_prcomp(MTBLS579, center = TRUE)
print_mpca(MTBLS_mpca)
}
