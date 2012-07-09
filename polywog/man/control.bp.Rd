\name{control.bp}
\alias{control.bp}
\title{Auxiliary for bootstrap options}
\usage{
  control.bp(reuse.lambda = FALSE, reuse.penwt = FALSE,
    maxtries = 1000, report = FALSE, scad.maxit = 5000)
}
\arguments{
  \item{reuse.lambda}{logical: whether to use the
  penalization parameter from the original fit
  (\code{TRUE}), or to cross-validate within each iteration
  (\code{FALSE}, default).}

  \item{reuse.penwt}{logical: whether to use the penalty
  weights from the original dataset for adaptive LASSO
  models (\code{TRUE}), or to re-calculate penalty weights
  within each iteration (\code{FALSE}, default).}

  \item{maxtries}{maximum number of attempts to generate a
  bootstrap sample with a non-collinear model matrix (often
  problematic with lopsided binary regressors) before
  failing.}

  \item{report}{logical: whether to print a status bar.
  Not available if \code{.parallel = TRUE}.}

  \item{scad.maxit}{maximum number of iterations for
  \code{\link{ncvreg}} in SCAD models.}
}
\value{
  A list containing the function arguments.
}
\description{
  Used to pass options to \code{\link{bootPolywog}} when
  bootstrapping via the \code{boot} option of the main
  \code{\link{polywog}} function.
}
\author{
  Brenton Kenkel and Curtis S. Signorino
}

