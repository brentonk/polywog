# polywog package news

## polywog 0.4-2 (May 2025)

* Under-the-hood changes to pass CRAN checks


## polywog 0.4-1 (April 2018)

* Dependency changes:
    * Occupational prestige data for the examples is now drawn from `carData` package (formerly `car`)

* Other minor changes to namespace calls to pass CRAN checks


## polywog 0.4-0 (April 2014)

* `predVals()` has been rewritten to compute fitted values according to the
  "observed value" approach advocated in the following article:

  > Hanmer, M. J. and Ozan Kalkan, K. (2013), Behind the Curve: Clarifying the
  > Best Approach to Calculating Predicted Probabilities and Marginal Effects
  > from Limited Dependent Variable Models. *American Journal of Political
  > Science*, 57: 263--277. doi: 10.1111/j.1540-5907.2012.00602.x

  For more details, see `?predVals`.

* k-fold cross-validation can now be performed in parallel for adaptive LASSO
  models, controlled via the `.parallel` argument of `polywog()`.  Parallel
  computation of bootstrap iterations is now handled via the `.parallel`
  argument of `bootPolywog()` or `control.bp()`.

* Adds `model.matrix` and `model.frame` methods for objects of class
  `"polywog"`

* Polynomial expansions of the design matrix are now handled in C++, and
  obsolete functions `polym2()` and `rawpoly()` have been removed

* New arguments of `polywog()`:
    * `lambda`, `nlambda`, and `lambda.min.ratio` for finer control of the
      sequence of penalization factor values examined
    * `foldid` for direct specification of cross-validation folds (only
      available when fitting via the adaptive LASSO)
    * `thresh` and `maxit` for finer control of the convergence criterion,
      replacing old argument `scad.maxit`

* Dependency changes:
    * All dependencies imported instead of attached to the search path, except
      `miscTools` which must be attached to provide the `margEff` generic
    * `glmnet` 1.9-5 required (for parallel cross-validation)
    * `ncvreg` 2.4-0 required (for bug fix in `cv.ncvreg`)
    * `iterators` and `Rcpp` required
    * `car` no longer required (but still suggested)
    * `matrixStats` and `games` no longer required


## polywog 0.3-0 (January 2013)

* `polywog()` now has argument `unpenalized` to exclude some terms from the
  adaptive LASSO penalty

* `bootPolywog()` now has argument `maxtries` to control failure when a
  non-collinear bootstrap model matrix cannot be found

* `bootPolywog()` now has argument `min.prop` to ensure a minimum amount of
  variation in the bootstrapped response variable in binary models

* The `fitted.values` element of `"polywog"` objects is now on the response
  scale instead of the link scale (i.e., transformed to probabilities when
  `family = "binomial"`)

* Fixed bug where the `polywog.fit` element of `cv.polywog()` output would not
  contain fitted values

* Fixed bug that sometimes caused `predVals()` to fail unexpectedly


## polywog 0.2-0 (June 2012)

* New function `cv.polywog()` to select both the polynomial degree and the
  penalization parameter by cross-validation

* New method `margEff.polywog()` to compute observation-wise and average
  marginal effects from a fitted model

* `varNames` element of a `"polywog"` object is now a character vector rather
  than a list (and is generated more safely)

* `"polyTerms"` attribute of matrix returned by `polym2()` is now a matrix
  rather than a data frame

* `predict.polywog()` now works correctly when `newdata` is a model frame


## polywog 0.1-0 (May 2012)

* Initial release
