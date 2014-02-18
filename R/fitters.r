##
## Calculate adaptive lasso results.  k-fold cross-validation is used to select
## the penalization parameter if 'lambda' is NULL.
##
## This function is used within 'polywog' (and 'bootPolywog'); it typically
## should not be called directly by a user.
##
## ARGUMENTS:
##   X: polynomial-expanded covariate matrix (without intercept)
##   y: response variable
##   weights: vector of observation weights
##   family: "gaussian" or "binomial"
##   penwt: vector of covariate weights for adaptive penalization
##   lambda: penalization parameter (NULL to cross-validate)
##   nfolds: number of cross-validation folds
##
## RETURN:
##   coef: coefficients from the fitted adaptive lasso model
##   lambda: penalization parameter used to obtain coefficients
##
fitALasso <- function(X, y, weights, family, penwt, lambda, nfolds, ...)
{
    if (is.null(lambda)) {
        ## Cross-validate
        res <- cvALasso(X = X, y = y, weights = weights, family = family,
                        penwt = penwt, nfolds = nfolds)
        coef <- res$coef
        lambda <- res$lambda
    } else {
        ## Run 'glmnet' directly
        res <- glmnet(x = X, y = y, weights = weights, family = family,
                      standardize = FALSE, penalty.factor = penwt,
                      lambda = lambda)
        coef <- c(res$a0[1], res$beta[, 1])
        names(coef)[1] <- "(Intercept)"
    }
    return(list(coef = coef, lambda = lambda))
}

##
## Analogue of fitALasso, for SCAD
##
fitSCAD <- function(X, y, family, lambda, nfolds, scad.maxit, ...)
{
    if (is.null(lambda)) {
        ## Cross-validate
        res <- cvSCAD(X = X, y = y, family = family, nfolds = nfolds,
                      scad.maxit = scad.maxit)
        coef <- res$coef
        lambda <- res$lambda
    } else {
        ## Run 'ncvreg' directly
        res <- ncvreg(X, y, family = family, penalty = "SCAD",
                      lambda = c(0, lambda), max.iter = scad.maxit)
        coef <- res$beta[, 2]
    }
    return(list(coef = coef, lambda = lambda))
}

cvALasso <- function(X, y, weights, family, penwt, nfolds, foldid = NULL, ...)
{
    ## Calling argument from a list to deal with NULL 'foldid' (the relevant
    ## flag in 'cv.glmnet' only deals with missing 'foldid')
    arglist <- list(x = X, y = y, weights = weights, family = family,
                    standardize = FALSE, penalty.factor = penwt,
                    nfolds = nfolds)
    if (!is.null(foldid))
        arglist$foldid <- foldid
    ans.cv <- do.call("cv.glmnet", arglist)

    ## Extract relevant info from fitted cross-validation object
    lambda <- ans.cv$lambda.min
    which.lambda <- which(ans.cv$lambda == lambda)
    coef <- c(ans.cv$glmnet.fit$a0[which.lambda],
              ans.cv$glmnet.fit$beta[, which.lambda])
    names(coef)[1] <- "(Intercept)"
    list(coef = coef, lambda = lambda, cverr = ans.cv$cvm[which.lambda])
}

cvSCAD <- function(X, y, family, nfolds, foldid = NULL, scad.maxit, ...)
{
    ans.cv <- ..cv.ncvreg(X = X, y = y, family = family, penalty = "SCAD",
                          max.iter = scad.maxit, nfolds = nfolds,
                          foldid = foldid)
    lambda <- ans.cv$lambda[ans.cv$min]
    which.lambda <- ans.cv$min
    ans <- ncvreg(X = X, y = y, family = family, penalty = "SCAD",
                  max.iter = scad.maxit)
    coef <- ans$beta[, which.lambda]
    list(coef = coef, lambda = lambda, cverr = ans.cv$cve[which.lambda])
}

##
## The following functions are taken from 'ncvreg' v2.3.2, used under GPL
## license.  'cv.ncvreg' is modified to allow user-specified fold IDs and to fix
## a bug in randomization for the Gaussian case.  The others are functions not
## exported from the 'ncvreg' namespace.
##

..cv.ncvreg <- function (X, y, family = c("gaussian", "binomial"), alpha = 1, 
    lambda.min = ifelse(n > p, 0.001, 0.05), nlambda = 100, lambda, 
    nfolds = 10, foldid, seed, trace = FALSE, ...) 
{
    family <- match.arg(family)
    if (alpha <= 0) 
        stop("alpha must be greater than 0; choose a small positive number instead")
    if (!missing(seed)) 
        set.seed(seed)
    n <- length(y)
    meanx <- apply(X, 2, mean)
    normx <- sqrt(apply((t(X) - meanx)^2, 1, sum)/n)
    nz <- which(normx > 1e-04)
    XX <- scale(X[, nz], meanx[nz], normx[nz])
    p <- ncol(XX)
    if (family == "gaussian") 
        yy <- y - mean(y)
    else yy <- y
    if (missing(lambda)) {
        lambda <- ..setupLambda(XX, yy, family, alpha, lambda.min, 
            nlambda)
        user.lambda <- FALSE
    }
    else {
        nlambda <- length(lambda)
        user.lambda <- TRUE
    }
    rm(XX)
    error <- array(NA, dim = c(nfolds, length(lambda)))

    ## (Not sure why two different algorithms originally, maybe just to ensure
    ## that there was variation in the outcome in all folds?  Either way, just
    ## going to use the same code as in 'cv.glmnet')
    ## 
    ## if (family == "gaussian") {
    ##     ## cv.ind <- ceiling((1:n)/n * nfolds)
    ##     cv.ind <- ceiling(sample(1:n)/n * nfolds)
    ## }
    ## else if (family == "binomial") {
    ##     ind1 <- which(y == 1)
    ##     ind0 <- which(y == 0)
    ##     n1 <- length(ind1)
    ##     n0 <- length(ind0)
    ##     cv.ind1 <- ceiling(sample(1:n1)/n1 * nfolds)
    ##     cv.ind0 <- ceiling(sample(1:n0)/n0 * nfolds)
    ##     cv.ind <- numeric(n)
    ##     cv.ind[y == 1] <- cv.ind1
    ##     cv.ind[y == 0] <- cv.ind0
    ## }
    if (missing(foldid) || is.null(foldid)) {
        cv.ind <- sample(rep(seq(nfolds), length = n))
    } else {
        cv.ind <- foldid
        nfolds <- max(foldid)
    }

    for (i in 1:nfolds) {
        if (trace) 
            cat("Starting CV fold #", i, sep = "", "\n")
        X1 <- X[cv.ind != i, ]
        y1 <- y[cv.ind != i]
        X2 <- X[cv.ind == i, ]
        y2 <- y[cv.ind == i]
        fit.i <- ncvreg(X1, y1, family = family, alpha = alpha, 
            lambda = lambda, warn = FALSE, ...)
        yhat <- predict(fit.i, X2, type = "response")
        error[i, 1:ncol(yhat)] <- ..loss.ncvreg(y2, yhat, family)
    }
    ind <- which(apply(is.finite(error), 2, all))
    E <- error[, ind]
    lambda <- lambda[ind]
    val <- list(E = E, cve = apply(E, 2, mean), lambda = lambda)
    val$min <- which.min(val$cve)
    class(val) <- "cv.ncvreg"
    return(val)
}

..setupLambda <- function (X, y, family, alpha, lambda.min, n.lambda) 
{
    n <- nrow(X)
    p <- ncol(X)
    if (family == "gaussian") {
        r <- y - mean(y)
        l1.max <- max(abs(crossprod(X, r)/n))
    }
    if (family == "binomial") {
        fit <- glm(y ~ 1, family = "binomial")
        pi. <- fit$fitted.values
        w <- pi. * (1 - pi.)
        r = (y - pi.)/w
        l1.max <- max(abs(crossprod(X, w * r)/n))
    }
    lambda.max <- l1.max/alpha
    if (lambda.min == 0) 
        lambda <- c(exp(seq(log(lambda.max), log(0.001 * lambda.max), 
            len = n.lambda - 1)), 0)
    else lambda <- exp(seq(log(lambda.max), log(lambda.min * 
        lambda.max), len = n.lambda))
    return(lambda)
}

..loss.ncvreg <- function (y, yhat, family) 
{
    n <- length(y)
    if (family == "gaussian") 
        return(apply(y - yhat, 2, crossprod)/(2 * length(y)))
    if (family == "binomial") 
        return((apply(log(yhat[y == 1, , drop = FALSE]), 2, sum) + 
            apply(log(1 - yhat[y == 0, , drop = FALSE]), 2, sum))/(-1 * 
            length(y)))
}

