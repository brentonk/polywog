\name{cv.polywog}
\alias{cv.polywog}
\title{Cross-validation for polywog models}
\usage{
  cv.polywog(formula, ..., method = c("alasso", "scad"),
    model = TRUE, X = FALSE, y = FALSE, degrees.cv = 1:3,
    nfolds = 10, scad.maxit = 5000)
}
\arguments{
  \item{formula}{model formula specifying the response and
  input variables.}

  \item{...}{other arguments to be passed to
  \code{\link{polywog}}.  Arguments controlling the
  bootstrap will be ignored.}

  \item{method}{variable selection method: \code{"alasso"}
  (default) for adaptive LASSO or \code{"scad"} for SCAD.}

  \item{model}{logical: whether to include the model frame
  in the returned object (default \code{TRUE}).  It may be
  problematic to run \code{\link{predVals}} or
  \code{\link{bootPolywog}} on a polywog object fit with
  \code{model = FALSE}.}

  \item{X}{logical: whether to include the
  polynomial-expanded design matrix in the returned object
  (default \code{FALSE}).}

  \item{y}{logical: whether to include the response
  variable in the returned object (default \code{FALSE}).}

  \item{degrees.cv}{vector of polynomial degrees to examine
  via cross-validation.}

  \item{nfolds}{number of folds to use in cross-validation
  to select the penalization parameter.}

  \item{scad.maxit}{maximum number of iterations when
  \code{method = "scad"} (see \code{\link{ncvreg}});
  ignored when \code{method = "alasso"}.}
}
\value{
  An object of class \code{"cv.polywog"}, a list
  containing: \describe{ \item{\code{results}}{A table of
  each degree tested, the optimal penalization parameter
  \eqn{\lambda} for that degree, and its cross-validation
  error.} \item{\code{degree.min}}{The polynomial degree
  giving the lowest cross-validation error.}
  \item{\code{polywog.fit}}{A \code{\link{polywog}} model,
  fit at the polynomial degree giving the lowest
  cross-validation error.} }

  Because the returned object contains the fitted polywog
  model for the optimal degree, no additional runs of
  \code{\link{polywog}} are necessary to estimate
  coefficients or the penalization parameter \eqn{\lambda}.
  However, bootstrap coefficients must be obtained by
  running \code{\link{bootPolywog}} on the
  \code{"polywog.fit"} element of the returned object, as
  in the examples below.
}
\description{
  Uses k-fold cross-validation to select the polynomial
  degree and penalization parameter for a
  \code{\link{polywog}} model.
}
\examples{
## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Examine degrees 1 through 4
set.seed(39)
cv1 <- cv.polywog(prestige ~ education + income + type, data = Prestige,
                  degrees.cv = 1:4, nfolds = 10)

print(cv1)

## Extract best model and bootstrap
fit1 <- cv1$polywog.fit
fit1 <- bootPolywog(fit1, nboot = 10)
summary(fit1)
}
\author{
  Brenton Kenkel and Curtis S. Signorino
}

