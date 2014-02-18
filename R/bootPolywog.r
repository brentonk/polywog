##' Auxiliary for bootstrap options
##'
##' Used to pass options to \code{\link{bootPolywog}} when bootstrapping via the
##' \code{boot} option of the main \code{\link{polywog}} function.
##' @inheritParams bootPolywog
##' @return A list containing the function arguments.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
control.bp <- function(reuse.lambda = FALSE, reuse.penwt = FALSE,
                       maxtries = 1000, min.prop = 0, report = FALSE,
                       scad.maxit = 5000)
{
    list(reuse.lambda = reuse.lambda, reuse.penwt = reuse.penwt,
         maxtries = maxtries, report = report, scad.maxit = scad.maxit)
}

##
## Innards of the bootstrap procedure
##
bootFit <- function(X, y, weights, family, lambda, penwt, method, penwt.method,
                    unpenalized, nfolds, scad.maxit, pb, i)
{
    ## Calculate penalty weights unless specified
    if (method == "alasso" && is.null(penwt)) {
        if (family == "gaussian" || penwt.method == "lm") {
            penwt <- lsfit(X, y, wt = weights, intercept = TRUE)$coef
            penwt <- 1 / abs(penwt[-1])
        } else {
            penwt <- penaltyWeightsBinary(X, y, weights)
        }

        penwt[unpenalized] <- 0
    }

    ## Fit model on bootstrap data and catch any convergence errors
    fitPolywog <- switch(method, alasso = fitALasso, scad = fitSCAD)
    ans <- tryCatch(fitPolywog(X = X, y = y, weights = weights, family = family,
                               penwt = penwt, lambda = lambda, nfolds = nfolds,
                               scad.maxit = scad.maxit)$coef,
                    error = identity)

    ## Increment progress bar if specified
    if(!is.null(pb))
        setTxtProgressBar(pb, i)

    return(ans)
}

##' Bootstrap a fitted polywog model
##'
##' Nonparametric bootstrap of the \code{\link{polywog}} regression procedure.
##' Can be run on a fitted model of class \code{"polywog"}, or within the
##' original procedure via the \code{boot} argument.
##'
##' When \code{.parallel = TRUE}, parallel computation is performed via
##' \code{\link[foreach:foreach]{\%dopar\%}} using the currently registered
##' backend.  Typically this will be \pkg{doMC} on Mac/Unix, \pkg{doSMP} on
##' Windows, and \pkg{doSNOW} in cluster environments.  Users must load the
##' appropriate packages and register the parallel environment before calling
##' \code{bootPolywog} (or \code{\link{polywog}} with \code{boot > 0}).  If a
##' parallel backend is not registered but \code{.parallel = TRUE}, computation
##' will proceed sequentially and \code{\%dopar\%} will issue a warning.
##' @param model a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param nboot number of bootstrap iterations.
##' @param reuse.lambda logical: whether to use the penalization parameter from
##' the original fit (\code{TRUE}), or to cross-validate within each iteration
##' (\code{FALSE}, default).
##' @param reuse.penwt logical: whether to use the penalty weights from the
##' original dataset for adaptive LASSO models (\code{TRUE}), or to re-calculate
##' penalty weights within each iteration (\code{FALSE}, default).
##' @param maxtries maximum number of attempts to generate a bootstrap sample
##' with a non-collinear model matrix (often problematic with lopsided binary
##' regressors) before failing.
##' @param min.prop for models with a binary response, minimum proportion of
##' non-modal outcome to ensure is included in each bootstrap iteration (for
##' example, set \code{min.prop = 0.1} to throw out any bootstrap iteration
##' where less than 10 percent of the responses are 1's)
##' @param report logical: whether to print a status bar.  Not available if
##' \code{.parallel = TRUE}.
##' @param scad.maxit maximum number of iterations for \code{\link{ncvreg}} in
##' SCAD models.
##' @param .parallel logical: whether to parallelize computation via
##' \code{\link[foreach]{foreach}}; see "Details" below.
##' @param .matrixOnly logical: whether to return just the matrix of bootstrap
##' coefficients (\code{TRUE}), or the originally supplied model with the
##' bootstrap matrix as the \code{boot.matrix} element (\code{FALSE}, default).
##' @return If \code{.matrixOnly = FALSE}, the returned object is \code{model}
##' with the bootstrap matrix included as its \code{boot.matrix} element.  If
##' \code{.matrixOnly = TRUE}, just the matrix is returned.  In either case, the
##' bootstrap matrix is a sparse matrix of class
##' \code{"\link[=dgCMatrix-class]{dgCMatrix}"}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @examples
##' ## Using occupational prestige data
##' data(Prestige, package = "car")
##' Prestige <- transform(Prestige, income = income / 1000)
##'
##' ## Fit a polywog model without bootstrap iterations
##' fit1 <- polywog(prestige ~ education + income + type, data = Prestige)
##' summary(fit1)
##'
##' ## Add bootstrap to the fitted model
##' fit2 <- bootPolywog(fit1, nboot = 10)
##' summary(fit2)
##'
##' ## Example of parallel processing on Mac/Unix via 'doMC'
##' \dontrun{
##' library(doMC)
##' registerDoMC()
##' 
##' fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)
##' }
##'
##' ## Example of parallel processing on Windows via 'doSMP'
##' \dontrun{
##' library(doSMP)
##' w <- startWorkers()
##' registerDoSMP(w)
##' 
##' fit2 <- bootPolywog(fit1, nboot = 100, .parallel = TRUE)
##' 
##' stopWorkers(w)
##' }
bootPolywog <- function(model, nboot = 100, reuse.lambda = FALSE,
                        reuse.penwt = FALSE, maxtries = 1000, min.prop = 0,
                        report = FALSE,
                        scad.maxit = 5000,
                        .parallel = FALSE, .matrixOnly = FALSE)
{
    ## Can't have a progress bar in parallel, sadly
    if (.parallel && report) {
        report <- FALSE
        warning("'report' set to FALSE because parallelization enabled")
    }

    ## Extract relevant information about initial fit
    ncf <- length(coef(model))
    formula <- model$formula
    degree <- model$degree
    family <- model$family
    method <- model$method
    penwt.method <- model$penwt.method
    unpenalized <- model$unpenalized
    pivot <- model$pivot
    nobs <- model$nobs
    lambda <- if (reuse.lambda) model$lambda else NULL
    penwt <- if (reuse.penwt) model$penwt else NULL
    nfolds <- model$nfolds
    weights <- model$weights
    if (is.null(weights))
        weights <- rep(1, nobs)

    ## Obtain original model matrix and response
    if (!is.null(model$X)) {
        X <- model$X
    } else {
        if (is.null(model$model))
            stop("Fitted object must contain either 'model' or both 'X' and 'y'; re-run polywog with \"model = TRUE\"")
        X <- makeX(formula, model$model, degree)[, pivot, drop = FALSE]
    }
    if (!is.null(model$y)) {
        y <- model$y
    } else {
        if (is.null(model$model))
            stop("Fitted object must contain either 'model' or both 'X' and 'y'; re-run polywog with \"model = TRUE\"")
        y <- model.part(formula, model$model, lhs = 1, drop = TRUE)
    }
    isBinary <- length(unique(y)) <= 2

    ## Bootstrap iterations
    pb <- if (report) txtProgressBar(min = 0, max = nboot) else NULL

    ## Loop over bootstrap iterations
    if (.parallel) {
        ## Loop in parallel via 'foreach'
        ans <- foreach (i = seq_len(nboot), .packages = "polywog") %dopar% {
            tries <- 0
            repeat {
                ## Ensure that the bootstrap X matrix has the same rank as the
                ## original one
                tries <- tries + 1
                if (tries > maxtries)
                    stop("'maxtries' reached; no non-collinear bootstrap sample found")
                ind <- sample(seq_len(nobs), nobs, replace = TRUE)
                Xgood <- qr(cbind(1L, X[ind, , drop = FALSE]))$rank == ncf
                ygood <- !isBinary || (mean(y[ind]) > min.prop &&
                                       mean(1-y[ind]) > min.prop)
                if (Xgood && ygood)
                    break
            }
            polywog:::bootFit(X = X[ind, , drop = FALSE], y = y[ind], weights =
                    weights[ind], family = family, lambda = lambda, penwt =
                    penwt, method = method, penwt.method = penwt.method,
                    unpenalized = unpenalized,
                    nfolds = nfolds, scad.maxit = scad.maxit, pb = pb, i = i)
        }
    } else {
        ## Use a sequential 'for' loop
        ans <- vector("list", nboot)
        for (i in seq_len(nboot)) {
            tries <- 0
            repeat {
                ## Ensure that the bootstrap X matrix has the same rank as the
                ## original one
                tries <- tries + 1
                if (tries > maxtries)
                    stop("'maxtries' reached; no non-collinear bootstrap sample found")
                ind <- sample(seq_len(nobs), nobs, replace = TRUE)
                if (qr(cbind(1L, X[ind, , drop = FALSE]))$rank == ncf)
                    break
            }
            ans[[i]] <-
                bootFit(X = X[ind, , drop = FALSE], y = y[ind], weights =
                        weights[ind], family = family, lambda = lambda, penwt =
                        penwt, method = method, penwt.method = penwt.method,
                        unpenalized = unpenalized,
                        nfolds = nfolds, scad.maxit = scad.maxit, pb = pb, i = i)
        }
    }

    if (report)  # Print newline after progress bar completes
        cat("\n")

    ## Warn about failures
    failures <- sapply(ans, function(x) inherits(x, "error"))
    if (sum(failures) > 0) {
        warning(sum(failures),
                " bootstrap iterations failed to converge and were removed",
                if (reuse.lambda) "; try setting 'reuse.lambda' to FALSE")
    }

    ## Store results in sparse matrix
    ans <- do.call(rbind, ans[!failures])
    ans <- Matrix(ans, sparse = TRUE)

    ## If .matrixOnly (typically only used within 'polywog'), return only the
    ## bootstrap matrix itself; otherwise, return the original model object with
    ## the bootstrap matrix included as an element
    if (.matrixOnly) {
        return(ans)
    } else {
        model$boot.matrix <- ans
        return(model)
    }
}
