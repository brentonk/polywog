## Undocumented "helper" functions

##
## Calculate penalty weights via logistic regression for adaptive lasso with
## binary outcomes
##
## ARGUMENTS:
##   X: polynomial expansion of covariate matrix (without intercept)
##   y: response variable
##   weights: vector of observation weights
##
## RETURN:
##   numeric vector of penalty weights for coefficients
##
penaltyWeightsBinary <- function(X, y, weights)
{
    ## Compute coefficients
    ans <- suppressWarnings(glm.fit(x = cbind(1L, X), y = y, weights = weights,
                                    family = binomial()))

    ## Separation check and convergence check (same as in 'glm.fit', but with
    ## warning messages tailored for this use case)
    if (!ans$converged)
        warning("'glm.fit' did not converge when computing penalty weights; consider using penwt.method = \"lm\"")
    eps <- 10 * .Machine$double.eps
    if (any(ans$fitted > 1 - eps) || any(ans$fitted < eps))
        warning("fitted probabilities numerically 0 or 1 occurred when computing penalty weights; consider using penwt.method = \"lm\"")

    ## Calculate weights from coefficients (exclude intercept)
    ans <- 1 / abs(ans$coef[-1])
    return(ans)
}
