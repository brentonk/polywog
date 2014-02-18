\name{predict.polywog}
\alias{predict.polywog}
\title{Predict method for polywog objects}
\usage{
\method{predict}{polywog}(object, newdata, type = c("link", "response"),
  interval = FALSE, level = 0.95, bag = FALSE, na.action = na.pass, ...)
}
\arguments{
  \item{object}{a fitted model of class \code{"polywog"},
  typically the output of \code{\link{polywog}}.}

  \item{newdata}{an optional data frame containing
  observations for which fitted values should be computed.
  If not specified, fitted values are generated for the
  data used to fit the model.}

  \item{type}{specifies whether the fitted values should be
  generated on the link scale (\eqn{X \beta}) or in terms
  of the expected value of the response variable.  These
  only differ for binomial family models.}

  \item{interval}{logical: whether to calculate bootstrap
  confidence intervals for each fitted value.}

  \item{level}{confidence level for the intervals.}

  \item{bag}{logical: whether to use "bootstrap
  aggregation" to generate the main fitted values (if
  \code{FALSE}, they are calculated from the main model
  fit).}

  \item{na.action}{a function specifying what to do with
  observations in \code{newdata} containing \code{NA}s
  (default \code{\link{na.pass}}).  See "Details".}

  \item{...}{other arguments, currently ignored.}
}
\value{
If \code{interval = TRUE}, a matrix containing each fitted
value and its confidence interval.  Otherwise, a vector
containing the fitted values.
}
\description{
Generates fitted values, including bootstrap confidence
intervals, for in- and out-of-sample data from a fitted
polywog model.
}
\details{
There are some special considerations involving elements of
\code{newdata} that contain \code{NA}s.  If all estimated
coefficients involving the variable \code{z} are exactly 0
(i.e., \code{z} is completely excluded from the estimated
model), then it is still possible to generate a fitted
value for an observation where \code{z} is not recorded.
Similarly, if all coefficients involving \code{z} are
estimated as 0 in every bootstrap iteration, it is still
possible to compute confidence intervals for such an
observation.

By default, \code{predict.polywog} computes fitted values
and confidence intervals whenever it is possible to do so.
To compute these values only for fully complete cases in
\code{newdata}, set \code{na.action = \link{na.omit}} or
\code{na.action = \link{na.exclude}}.  \code{na.omit} will
return an object whose dimension equals the number of
complete cases; \code{na.exclude} gives the same dimension
as \code{newdata}, padded with \code{NA}s where
appropriate.
}
\author{
Brenton Kenkel and Curtis S. Signorino
}
\seealso{
For more user-friendly generation of fitted values, see
\code{\link{predVals}}.  To compute marginal effects, see
\code{\link{margEff.polywog}}.
}

