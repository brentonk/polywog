\name{predVals}
\alias{predVals}
\title{Easy computation of fitted values}
\usage{
predVals(model, xvars, xlims = list(), n = 100, interval = TRUE,
  level = 0.95, bag = FALSE, maxrows = 10000, ...)
}
\arguments{
  \item{model}{a fitted model of class \code{"polywog"},
  typically the output of \code{\link{polywog}}.}

  \item{xvars}{a character vector containing names of raw
  input variables (from \code{model$varNames}).  Partial
  matches are allowed.}

  \item{xlims}{named list of limits for the evaluation grid
  for each continuous variable in \code{xvars}.  If not
  given, the variable's observed range is used.}

  \item{n}{number of grid points at which to evaluate each
  continuous variable in \code{xvars}.}

  \item{interval}{logical: whether to compute bootstrap
  confidence intervals for each fitted value.}

  \item{level}{confidence level for the intervals.}

  \item{bag}{logical: whether to use "bootstrap
  aggregation" to generate the main fitted values (if
  \code{FALSE}, they are calculated from the main model
  fit).}

  \item{maxrows}{maximum number of rows of output.  Used to
  prevent accidental memory overruns when \code{xvars}
  contains more than two continuous variables.}

  \item{...}{used to set values for the variables other
  than \code{xvars} in the profile of observations.  See
  "Details" below.}
}
\value{
A data frame containing the fitted values and confidence
intervals (if requested) for each combination of covariate
values.
}
\description{
User-friendly generation of fitted values and their
confidence intervals from models of class \code{"polywog"},
using \code{\link{predict.polywog}} as a backend.
}
\details{
\code{predVals} allows users to examine the estimated
effects of input variables on the expected outcome using
the coefficients returned by \code{\link{polywog}}.  The
procedure is designed so that, for a preliminary analysis,
the user can simply specify the fitted model and the
independent variable of interest, and quickly obtain
predicted values.  However, it is flexible enough to allow
for finely tuned analysis as well.  The function is very
similar to \code{\link[games]{predProbs}} in the
\pkg{games} package.

The procedure works by varying \code{xvars}, the variables
of interest, across their observed ranges (or those
specified by the user in \code{xlims}) while holding all
other independent variables in the model fixed.  The
profile created by default is as follows (the same defaults
as in the \code{\link[Zelig]{sim}} function in the
\pkg{Zelig} package): \itemize{ \item numeric, non-binary
variables are fixed at their means \item
\code{\link{ordered}} and binary variables are fixed at
their medians \item all others are fixed at their modes
(see \code{\link[games]{Mode}}) } However, it is possible
to override these defaults for any or all variables.  For
example, to set a variable named \code{polity} to its lower
quartile, call \code{predVals} with the argument
\code{polity = quantile(polity, 0.25)}.  To set a factor
variable to a particular level, provide the name of the
level as a character string (in quotes).  See the examples
below for illustrations of this functionality.

All confidence intervals are generated via the bootstrap.
If \code{model} does not have a \code{boot.matrix} element
(see \code{\link{bootPolywog}}), confidence intervals will
not be computed.
}
\examples{
## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model with bootstrap iterations
set.seed(22)
fit1 <- polywog(prestige ~ education + income + type, data = Prestige,
                boot = 10)

## Predicted prestige across occupational categories
predVals(fit1, "type")
predVals(fit1, "type", income = quantile(income, 0.25))

## Predicted prestige by education
predVals(fit1, "education", n = 10)
predVals(fit1, "education", n = 10, income = quantile(income, 0.25))
}
\author{
Brenton Kenkel and Curtis S. Signorino
}
\seealso{
\code{\link{predict.polywog}} for more flexible (but less
user-friendly) computation of fitted values.
\code{\link{plot.polywog}} for plotting fitted values and
their confidence intervals.
}

