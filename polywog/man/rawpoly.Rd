\name{rawpoly}
\alias{rawpoly}
\title{Calculate raw polynomials}
\usage{
  rawpoly(x, degree = 1, na.ok = FALSE)
}
\arguments{
  \item{x}{a numeric vector.}

  \item{degree}{the degree of the polynomial.}

  \item{na.ok}{logical: whether to allow \code{NA}s.  If
  \code{na.ok = FALSE}, then the function will stop with an
  error if any \code{NA}s are present in \code{x}.}
}
\value{
  A matrix containing a raw polynomial expansion of
  \code{x}.
}
\description{
  This is essentially the same as running \code{poly(x,
  degree = degree, raw = TRUE)} (and is based on the
  \code{\link{poly}} code), but \code{NA}s are allowed and
  no warning is issued when the degree exceeds the number
  of unique points.
}
\examples{
x <- c(1, 0, 1, 0, 0, NA, 1)
rawpoly(x, degree = 3, na.ok = TRUE)
}
\author{
  Brenton Kenkel and Curtis S. Signorino
}

