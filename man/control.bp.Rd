\name{control.bp}
\alias{control.bp}
\title{Auxiliary for bootstrap options}
\usage{
control.bp(.parallel = FALSE, reuse.lambda = FALSE, reuse.penwt = FALSE,
  nlambda = 100, lambda.min.ratio = 1e-04, nfolds = 10, thresh = NULL,
  maxit = NULL, maxtries = 1000, min.prop = 0, report = FALSE)
}
\arguments{
  \item{.parallel}{logical: whether to perform computations
  in parallel using a backend registered with
  \code{\link{foreach}}.}

  \item{reuse.lambda}{logical: whether to use the
  penalization parameter from the original fit
  (\code{TRUE}), or to cross-validate within each iteration
  (\code{FALSE}, default).}

  \item{reuse.penwt}{logical: whether to use the penalty
  weights from the original fit (\code{TRUE}), or to
  re-calculate them within each iteration (\code{FALSE},
  default).}

  \item{nlambda}{number of values of the penalty factor to
  examine in cross-validation, as in
  \code{\link{polywog}}.}

  \item{lambda.min.ratio}{ratio of the smallest value of
  the penalty factor to the largest, as in
  \code{\link{polywog}}.}

  \item{nfolds}{number of cross-validation folds to use.}

  \item{thresh}{convergence threshold, as in
  \code{\link{polywog}}.  If \code{NULL}, use the same
  value as in the original model.}

  \item{maxit}{iteration limit for fitting, as in
  \code{\link{polywog}}.  If \code{NULL}, use the same
  value as in the original model.}

  \item{maxtries}{maximum number of attempts to generate a
  bootstrap sample with a non-collinear model matrix (often
  problematic with lopsided binary regressors) before
  stopping and issuing an error message.}

  \item{min.prop}{for models with a binary response
  variable, minimum proportion of non-modal outcome to
  ensure is included in each bootstrap iteration (for
  example, set \code{min.prop = 0.1} to throw out any
  bootstrap iteration where less than 10 percent or more
  than 90 percent of the observations are 1's).}

  \item{report}{logical: whether to print a status bar.
  Not available if \code{.parallel = TRUE}.}
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

