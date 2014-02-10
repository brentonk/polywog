\name{polym2}
\alias{polym2}
\title{Calculate polynomial basis expansion of a matrix}
\usage{
polym2(x, degree = 1, na.ok = FALSE)
}
\arguments{
  \item{x}{a numeric matrix.}

  \item{degree}{the degree of the polynomial.}

  \item{na.ok}{logical: whether to allow \code{NA}s.  If
  \code{na.ok = FALSE}, then the function will stop with an
  error if any \code{NA}s are present in \code{x}.}
}
\value{
A numeric matrix containing the specified basis expansion
of the input matrix, along with attributes \code{degree}
giving the degree of each column and \code{polyTerms}
giving the power of each input term represented in each
column.
}
\description{
This produces the same output as running \code{polym(x,
degree = degree, raw = TRUE)} (and is based on the
\code{\link{polym}} code).  The differences are: \itemize{
\item Computation is much faster than \code{polym} when
\code{degree} is relatively low and \code{x} contains many
columns (the typical use case for \code{\link{polywog}}).
\item Column names of the output use the column names of
\code{x} in order to be more intelligible.  \item
\code{NA}s are allowed via the \code{na.ok} option. }
}
\examples{
x <- matrix(rnorm(20), ncol = 2)
colnames(x) <- c("a", "b")
polym2(x, degree = 3)
}
\author{
Brenton Kenkel
}

