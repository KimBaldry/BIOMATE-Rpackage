% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isempty.R
\name{is.empty}
\alias{is.empty}
\title{Check if an object is empty}
\usage{
is.empty(x, first.only = TRUE, all.na.empty = TRUE)
}
\arguments{
\item{x}{an object}

\item{first.only}{Are we only investigating the fist entry of a character vector, or all entries of a character vector (in this case a vector is returned)?}

\item{all.na.empty}{is a vector with all NA values counted as empty?}
}
\value{
Logical result or vector of logical
}
\description{
This is a function that checks if an object (charactor, vector, integer) has no information contained with it (including null, empty, length of 0).
}
\author{
Kimberlee Baldry
}
