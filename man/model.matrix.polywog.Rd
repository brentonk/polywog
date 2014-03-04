\name{model.matrix.polywog}
\alias{model.matrix.polywog}
\title{Model matrix of a polywog model}
\usage{
\method{model.matrix}{polywog}(object, type = c("raw", "expanded"), ...)
}
\arguments{
  \item{object}{a fitted model of class \code{"polywog"}}

  \item{type}{\code{"raw"}, the default, returns the
  non-expanded model matrix with no intercept (same number
  of columns as \code{object$polyTerms}).
  \code{"expanded"} returns the polynomial expansion used
  in fitting (number of columns equals
  \code{length(object$coefficients)}).}

  \item{...}{other arguments to be passed to further
  methods (typically only used internally)}
}
\value{
The design matrix of the specified model, consisting either
of raw terms or the full polynomial expansion depending on
the \code{type} argument.
}
\description{
Constructs the design matrix used to fit a
\code{\link{polywog}} model, similar to
\code{\link{model.matrix.lm}}.
}
\details{
There are two types of model matrix a user might want to
construct. First, there is the matrix of the raw input
terms that go into the eventual polynomial expansion.  Such
a matrix can be obtained by using \code{type = "raw"} (the
default).  To save memory, all of the fitting and analysis
functions in the \pkg{polywog} package use this form of the
design matrix, only constructing the polynomial expansion
on the fly as needed.

The other form of the model matrix is the full polynomial
expansion, where each column contains some power of the raw
inputs.  This can be obtained by using \code{type =
"expanded"}.  Depending on the number of observations, the
number of raw input variables, and the degree of the
polynomial expansion used in the fitted model, the
resulting matrix may require an unwieldy amount of memory.
}
\author{
Brenton Kenkel and Curtis S. Signorino
}

