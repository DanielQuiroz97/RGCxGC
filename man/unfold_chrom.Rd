% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unfold_chrom.R
\name{unfold_chrom}
\alias{unfold_chrom}
\title{Unfold two-dimensional chromatograms}
\usage{
unfold_chrom(Object)
}
\arguments{
\item{Object}{a batch_2DCOW or joined_chrom objects.}
}
\description{
`unfold` converts the two-dimensional chromatograms into a one dimensional
vector. Then, all chromatograms are joined into a matrix.
}
\details{
This function takes a single argument, batch_2DCOW or joined_chrom objects
and extracts each chromatogram and then it is unfolded into a one-dimensional
vector. Then, each one dimensional vector is joined in a single matrix, where
each row represent an observation or a chromatogram and each column
represent a variable, in our case, each retention time.
Also, in order to keep information about chromatographic runs, the retention
times for both dimensions are also exported.
}
\examples{
 

data(Myrothecium)
# Unfold 2D chromatogram
chrom_1D <- unfold_chrom(Myrothecium)
# Retrieve retention time for the first dimension
time_1D <- chrom_1D$time
# Retrieve the modulation time
modulation <- chrom_1D$mod_time
}
