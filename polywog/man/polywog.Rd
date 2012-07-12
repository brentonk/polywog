\name{polywog}
\alias{polywog}
\title{Polynomial regression with oracle variable selection}
\usage{
  polywog(formula, data, subset, weights, na.action,
    degree = 3, family = c("gaussian", "binomial"),
    method = c("alasso", "scad", "none"),
    penwt.method = c("lm", "glm"), penalize.linear = TRUE,
    boot = 0, control.boot = control.bp(),
    .parallel = FALSE, model = TRUE, X = FALSE, y = FALSE,
    nfolds = 10, scad.maxit = 5000)
}
\arguments{
  \item{formula}{model formula specifying the response and
  input variables.  See "Details" for more information.}

  \item{data}{a data frame, list or environment containing
  the variables specified in the model formula.}

  \item{subset}{an optional vector specifying a subset of
  observations to be used in fitting.}

  \item{weights}{an optional vector specifying weights for
  each observation to be used in fitting.}

  \item{na.action}{a function specifying what to do with
  observations containing \code{NA}s (default
  \code{\link{na.omit}}).}

  \item{degree}{integer specifying the degree of the
  polynomial expansion of the input variables.}

  \item{family}{\code{"gaussian"} (default) or
  \code{"binomial"} for logistic regression (binary
  response only).}

  \item{method}{variable selection method: \code{"alasso"}
  (default) for adaptive LASSO or \code{"scad"} for SCAD.
  You can also select \code{method = "none"} to return the
  model matrix and other information without fitting.}

  \item{penwt.method}{estimator for obtaining first-stage
  estimates in logistic models when \code{method =
  "alasso"}: \code{"lm"} (default) for a linear probability
  model, \code{"glm"} for logistic regression.}

  \item{penalize.linear}{logical: whether the penalty
  should be applied to any terms that are excluded from the
  basis expansion (see "Details" below).  Not available if
  \code{method = "scad"}.}

  \item{boot}{number of bootstrap iterations (0 for no
  bootstrapping).}

  \item{control.boot}{list of arguments to be passed to
  \code{\link{bootPolywog}} when bootstrapping; see
  \code{\link{control.bp}}.}

  \item{.parallel}{logical: whether to perform bootstrap
  iterations in parallel using
  \code{\link[foreach]{foreach}}.  See the "Details"
  section of the \code{\link{bootPolywog}} documentation
  page for more on parallel computation.}

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

  \item{nfolds}{number of folds to use in cross-validation
  to select the penalization parameter.}

  \item{scad.maxit}{maximum number of iterations when
  \code{method = "scad"} (see \code{\link{ncvreg}});
  ignored when \code{method = "alasso"}.}
}
\value{
  An object of class \code{"polywog"}, a list containing:
  \describe{ \item{\code{coefficients}}{the estimated
  coefficients.} \item{\code{lambda}}{the penalization
  parameter selected by cross-validation.}
  \item{\code{fitted.values}}{the fitted mean values for
  each observation used in fitting.}
  \item{\code{lmcoef}}{coefficients from an unpenalized
  least-squares regression.} \item{\code{penwt}}{penalty
  weights if \code{method = "alasso"}, \code{NULL}
  otherwise.} \item{\code{formula}}{model formula.}
  \item{\code{degree}}{degree of polynomial basis
  expansion.} \item{\code{family}}{model family,
  \code{"gaussian"} or \code{"binomial"}.}
  \item{\code{weights}}{observation weights if specified,
  \code{NULL} otherwise.}
  \item{\code{method}}{regularization method,
  \code{"alasso"} or \code{"scad"}.}
  \item{\code{penwt.method}}{estimator for penalty weights
  in adaptive LASSO models, \code{"lm"} or \code{"glm"}; is
  \code{NULL} if \code{method = "scad"}.}
  \item{\code{penalize.linear}}{whether terms included
  linearly are penalized.} \item{\code{which.linear}}{which
  columns of the model matrix (without intercept) are terms
  not included in the basis expansion.}
  \item{\code{nfolds}}{number of cross-validation folds.}
  \item{\code{terms}}{the \code{\link{terms}} object used
  in fitting.} \item{\code{pivot}}{indices of the
  non-collinear columns of the full basis expansion of the
  input variables.} \item{\code{nobs}}{number of
  observations used in fitting.}
  \item{\code{na.action}}{information on how \code{NA}s in
  the input data were handled.}
  \item{\code{xlevels}}{levels of factor variables used in
  fitting.} \item{\code{polyTerms}}{a matrix recording how
  many powers of each raw input term are represented in
  each column of the design matrix.}
  \item{\code{varNames}}{names of the raw input variables
  included in the model.} \item{\code{call}}{the original
  function call.} \item{\code{model}}{(if requested) the
  model frame (see \code{\link{model.frame}}).}
  \item{\code{X}}{(if requested) the basis-expanded design
  matrix.} \item{\code{y}}{(if requested) the response
  variable.} \item{\code{boot.matrix}}{(if \code{boot > 0})
  a sparse matrix of class
  \code{"\link[=dgCMatrix-class]{dgCMatrix}"} containing
  bootstrap coefficients (see \code{\link{bootPolywog}}).}
  }
}
\description{
  Fits a regression model using a polynomial basis
  expansion of the input variables, with penalization via
  the adaptive LASSO or SCAD to provide oracle variable
  selection.
}
\details{
  The design matrix for the regression is a polynomial
  basis expansion of the matrix of raw input variables.
  This includes all powers and interactions of the input
  variables up to the specified \code{degree}.  For
  example, the following terms will be included in
  \code{polywog(y ~ x1 + x2, degree = 3, ...)}: \itemize{
  \item terms of degree 0: intercept \item terms of degree
  1: \code{x1}, \code{x2} \item terms of degree 2:
  \code{x1^2}, \code{x2^2}, \code{x1*x2} \item terms of
  degree 3: \code{x1^3}, \code{x2^3}, \code{x1*x2^2},
  \code{x1^2*x2} } To exclude certain terms from the basis
  expansion, use a model formula like \code{y ~ x1 + x2 |
  z1 + z2}.  Only the degree 1 terms of \code{z1} and
  \code{z2} will be included.

  It is possible that the "raw" basis expansion will be
  rank-deficient, such as if there are binary input
  variables (in which case \eqn{x_i = x_i^n} for all \eqn{n
  > 0}).  The procedure detects collinearity via
  \code{\link{qr}} and removes extraneous columns before
  fitting.

  For both the adaptive LASSO and SCAD, the penalization
  factor \eqn{\lambda} is chosen by k-fold
  cross-validation.  The selected value minimizes the
  average mean squared error of out-of-sample fits.  (To
  select both \eqn{\lambda} and the polynomial degree
  simultaneously via cross-validation, see
  \code{\link{cv.polywog}}.)

  The bootstrap iterations may be run in parallel via
  \code{\link[foreach]{foreach}} by registering an
  appropriate backend and specifying \code{.parallel =
  TRUE}.  For more information on parallel processing, see
  the documentation and examples for
  \code{\link{bootPolywog}}.
}
\examples{
## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model with bootstrap iterations
set.seed(22)
fit1 <- polywog(prestige ~ education + income + type, data = Prestige,
                boot = 10)

## Basic information
print(fit1)
summary(fit1)

## See how fitted values change with education holding all else fixed
predVals(fit1, "education", n = 10)

## Plot univariate relationships
plot(fit1)

## Change regularization method
fit2 <- update(fit1, method = "scad")
cbind(coef(fit1), coef(fit2))
}
\author{
  Brenton Kenkel and Curtis S. Signorino
}
\references{
  Brenton Kenkel and Curtis S. Signorino.  2012.  "A Method
  for Flexible Functional Form Estimation: Bootstrapped
  Basis Regression with Variable Selection."  Typescript,
  University of Rochester.
}
\seealso{
  To estimate variation via the bootstrap, see
  \code{\link{bootPolywog}}.  To generate fitted values,
  see \code{\link{predVals}} (and the underlying method
  \code{\link{predict.polywog}}).  For plots, see
  \code{\link{plot.polywog}}. The polynomial degree may be
  selected via cross-validation using
  \code{\link{cv.polywog}}.

  Polynomial basis expansions of matrix inputs are computed
  with \code{\link{polym2}}.

  Adaptive LASSO estimates are provided via
  \code{\link{glmnet}} and \code{\link{cv.glmnet}} from the
  \pkg{glmnet} package.  SCAD estimates are via
  \code{\link{cv.ncvreg}} and \code{\link{ncvreg}} in the
  \pkg{ncvreg} package.
}

