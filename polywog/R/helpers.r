## Undocumented "helper" functions

##
## Remove all intercepts from a formula with a single left-hand side and one or
## more right-hand sides.
##
## ARGUMENTS:
##   formula.: an object of class "Formula"
##   mf: the model frame (included to deal with the case where 'formula.' is of
##       the form y ~ .)
##
## RETURN:
##   'formula.', with each rhs term including "- 1"
##
removeIntercepts <- function(formula., mf)
{
    nrhs <- length(formula.)[2]
    flist <- vector("list", nrhs)

    ## Keep LHS for first term [the first line is to ensure no error if the
    ## formula is of the form y ~ .]
    ft <- terms(formula(formula., lhs = 1, rhs = 1), data = mf)
    flist[[1]] <- update(formula(ft), . ~ . - 1)

    ## Don't use LHS in subsequent terms (if any)
    for (i in seq_len(nrhs-1)) {
        ft <- terms(formula(formula., lhs = 0, rhs = i+1), data = mf)
        flist[[i+1]] <- update(formula(ft), ~ . - 1)
    }

    ans <- do.call(as.Formula, flist)
    return(ans)
}

##
## Assemble polynomial-expanded model matrix (plus linear terms if specified)
## without intercept.
##
## This function is used within 'polywog' and some post-estimation analysis
## functions; it typically should not be called directly by a user.
##
## ARGUMENTS:
##   formula: model formula (of class "Formula")
##   mf: model frame
##   degree: degree of polynomial expansion
##   na.ok: whether to let NAs pass through
##
## RETURN:
##   model matrix
##
makeX <- function(formula, mf, degree, na.ok = FALSE)
{
    ## Ensure no intercepts included
    formula <- removeIntercepts(formula, mf)

    ## Polynomial expansion
    X <- model.matrix(formula, data = mf, rhs = 1)
    X <- polym2(X, degree = degree, na.ok = na.ok)

    ## Linear terms (if any)
    if (length(formula)[2] > 1) {
        ## Need to save "polyTerms" attribute or cbind() will destroy it
        pt1 <- attr(X, "polyTerms")
        Xlin <- model.matrix(formula, data = mf, rhs = 2)
        X <- cbind(X, Xlin)

        ## Add the new variables to "polyTerms"
        nlin <- ncol(Xlin)
        pt <- cbind(pt1, matrix(0, nrow = nrow(pt1), ncol = nlin))
        pt <- rbind(pt,
                    cbind(matrix(0, nrow = nlin, ncol = ncol(pt1)), diag(nlin)))
        colnames(pt) <- c(colnames(pt1), colnames(Xlin))
        attr(X, "polyTerms") <- pt
    } else {
        nlin <- 0
    }
    attr(X, "which.linear") <- rev(seq_len(ncol(X)))[seq_len(nlin)]

    return(X)
}

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
