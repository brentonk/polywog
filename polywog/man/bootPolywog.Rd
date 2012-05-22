\name{bootPolywog}
\alias{bootPolywog}
\title{Bootstrap a fitted polywog model}
\usage{
  bootPolywog(model, nboot = 100, reuse.lambda = FALSE,
    reuse.penwt = FALSE, report = FALSE, scad.maxit = 5000,
    .parallel = FALSE, .matrixOnly = FALSE)
}
\arguments{
  \item{model}{a fitted model of class \code{"polywog"},
  typically the output of \code{\link{polywog}}.}

  \item{nboot}{number of bootstrap iterations.}

  \item{reuse.lambda}{logical: whether to use the
  penalization parameter from the original fit
  (\code{TRUE}), or to cross-validate within each iteration
  (\code{FALSE}, default).}

  \item{reuse.penwt}{logical: whether to use the penalty
  weights from the original dataset for adaptive LASSO
  models (\code{TRUE}), or to re-calculate penalty weights
  within each iteration (\code{FALSE}, default).}

  \item{report}{logical: whether to print a status bar.
  Not available if \code{.parallel = TRUE}.}

  \item{scad.maxit}{maximum number of iterations for
  \code{\link{ncvreg}} in SCAD models.}

  \item{.parallel}{logical: whether to parallelize
  computation via \code{\link[foreach]{foreach}}; see
  "Details" below.}

  \item{.matrixOnly}{logical: whether to return just the
  matrix of bootstrap coefficients (\code{TRUE}), or the
  originally supplied model with the bootstrap matrix as
  the \code{boot.matrix} element (\code{FALSE}, default).}
}
\value{
  If \code{.matrixOnly = FALSE}, the returned object is
  \code{model} with the bootstrap matrix included as its
  \code{boot.matrix} element.  If \code{.matrixOnly =
  TRUE}, just the matrix is returned.  In either case, the
  bootstrap matrix is a sparse matrix of class
  \code{"\link[=dgCMatrix-class]{dgCMatrix}"}.
}
\description{
  Nonparametric bootstrap of the \code{\link{polywog}}
  regression procedure. Can be run on a fitted model of
  class \code{"polywog"}, or within the original procedure
  via the \code{boot} argument.
}
\details{
  When \code{.parallel = TRUE}, parallel computation is
  performed via \code{\link[foreach:foreach]{\%dopar\%}}
  using the currently registered backend.  Typically this
  will be \pkg{doMC} on Mac/Unix, \pkg{doSMP} on Windows,
  and \pkg{doSNOW} in cluster environments.  Users must
  load the appropriate packages and register the parallel
  environment before calling \code{bootPolywog} (or
  \code{\link{polywog}} with \code{boot > 0}).  If a
  parallel backend is not registered but \code{.parallel =
  TRUE}, computation will proceed sequentially and
  \code{\%dopar\%} will issue a warning.
}
\examples{
## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model without bootstrap iterations
fit1 <- polywog(prestige ~ education + income + type, data = Prestige)
summary(fit1)

## Add bootstrap to the fitted model
fit2 <- bootPolywog(fit1, nboot = 10)
summary(fit2)

## Example of parallel processing on Mac/Unix via 'doMC'
\dontrun{
library(doMC)
registerDoMC()

fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)
}

## Example of parallel processing on Windows via 'doSMP'
\dontrun{
library(doSMP)
w <- startWorkers()
registerDoSMP(w)

fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)

stopWorkers(w)
}
}
\author{
  Brenton Kenkel and Curtis S. Signorino
}

