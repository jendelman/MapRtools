% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LGtrim.R
\name{LGtrim}
\alias{LGtrim}
\title{Trim a linkage group based on genotype frequencies}
\usage{
LGtrim(geno, LODmat, thresh)
}
\arguments{
\item{geno}{matrix of haplotype dosages (markers x samples)}

\item{LODmat}{matrix of LOD scores for the markers}

\item{thresh}{numeric vector of thresholds for clusterings}
}
\value{
Either a ggplot2 object or a vector of marker names (see Details)
}
\description{
Trim a linkage group based on genotype frequencies
}
\details{
This function should only be run on a single linkage group (to form the linkage groups, use \code{\link{LG}}. If \code{thresh} is a numeric vector with multiple LOD thresholds, the function returns a plot showing the impact of the threshold on genotype frequencies. If \code{thresh} is a single value, the function returns a vector of the marker names that are retained. The rownames of \code{geno} and \code{LODmat} must match.
}
