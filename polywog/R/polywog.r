##' @include helpers.r
##' @include fitters.r
NULL

##' Bootstrapped basis regression with oracle model selection
##'
##' A package for flexible functional form estimation via bootstrapped basis
##' regression with oracle model selection.  This version of the software should
##' be considered \strong{in beta}.  Please email the maintainer, Brenton Kenkel
##' (\email{brenton.kenkel@@gmail.com}), with bug reports or feature requests.
##' @name polywog-package
##' @docType package
##' @section Acknowledgements: We thank the Wallis Institute of Political
##' Economy for financial support.
##' @references
##' Brenton Kenkel and Curtis S. Signorino.  2012.  "A Method for Flexible
##' Functional Form Estimation: Bootstrapped Basis Regression with Variable
##' Selection."  Typescript, University of Rochester.
NULL

##' Polynomial regression with oracle variable selection
##'
##' Fits a regression model using a polynomial basis expansion of the input
##' variables, with penalization via the adaptive LASSO or SCAD to provide
##' oracle variable selection.
##'
##' The design matrix for the regression is a polynomial basis expansion of the
##' matrix of raw input variables.  This includes all powers and interactions of
##' the input variables up to the specified \code{degree}.  For example, the
##' following terms will be included in \code{polywog(y ~ x1 + x2, degree = 3,
##' ...)}:
##' \itemize{
##'   \item terms of degree 0: intercept
##'   \item terms of degree 1: \code{x1}, \code{x2}
##'   \item terms of degree 2: \code{x1^2}, \code{x2^2}, \code{x1*x2}
##'   \item terms of degree 3: \code{x1^3}, \code{x2^3}, \code{x1*x2^2},
##' \code{x1^2*x2}
##' }
##' To exclude certain terms from the basis expansion, use a model formula like
##' \code{y ~ x1 + x2 | z1 + z2}.  Only the degree 1 terms of \code{z1} and
##' \code{z2} will be included.
##'
##' It is possible that the "raw" basis expansion will be rank-deficient, such
##' as if there are binary input variables (in which case \eqn{x_i = x_i^n} for
##' all \eqn{n > 0}).  The procedure detects collinearity via \code{\link{qr}} and
##' removes extraneous columns before fitting.
##'
##' For both the adaptive LASSO and SCAD, the penalization factor \eqn{\lambda}
##' is chosen by k-fold cross-validation.  The selected value minimizes the
##' average mean squared error of out-of-sample fits.  (To select both
##' \eqn{\lambda} and the polynomial degree simultaneously via cross-validation,
##' see \code{\link{cv.polywog}}.)
##'
##' The bootstrap iterations may be run in parallel via
##' \code{\link[foreach]{foreach}} by registering an appropriate backend and
##' specifying \code{.parallel = TRUE}.  For more information on parallel
##' processing, see the documentation and examples for \code{\link{bootPolywog}}.
##' @param formula model formula specifying the response and input
##' variables.  See "Details" for more information.
##' @param data a data frame, list or environment containing the variables
##' specified in the model formula.
##' @param subset an optional vector specifying a subset of observations to be
##' used in fitting.
##' @param weights an optional vector specifying weights for each observation to
##' be used in fitting.
##' @param na.action a function specifying what to do with observations
##' containing \code{NA}s (default \code{\link{na.omit}}).
##' @param degree integer specifying the degree of the polynomial expansion of
##' the input variables.
##' @param family \code{"gaussian"} (default) or \code{"binomial"} for logistic
##' regression (binary response only).
##' @param method variable selection method: \code{"alasso"} (default) for
##' adaptive LASSO or \code{"scad"} for SCAD.  You can also select \code{method
##' = "none"} to return the model matrix and other information without fitting.
##' @param penwt.method estimator for obtaining first-stage estimates in
##' logistic models when \code{method = "alasso"}: \code{"lm"} (default) for a
##' linear probability model, \code{"glm"} for logistic regression.
##' @param boot number of bootstrap iterations (0 for no bootstrapping).
##' @param control.boot list of arguments to be passed to
##' \code{\link{bootPolywog}} when bootstrapping; see \code{\link{control.bp}}.
##' @param .parallel logical: whether to perform bootstrap iterations in
##' parallel using \code{\link[foreach]{foreach}}.  See the "Details" section of
##' the \code{\link{bootPolywog}} documentation page for more on parallel
##' computation.
##' @param model logical: whether to include the model frame in the returned
##' object (default \code{TRUE}).  It may be problematic to run
##' \code{\link{predVals}} or \code{\link{bootPolywog}} on a polywog object fit
##' with \code{model = FALSE}.
##' @param X logical: whether to include the polynomial-expanded design matrix
##' in the returned object (default \code{FALSE}).
##' @param y logical: whether to include the response variable in the returned
##' object (default \code{FALSE}).
##' @param nfolds number of folds to use in cross-validation to select the
##' penalization parameter.
##' @param scad.maxit maximum number of iterations when \code{method = "scad"}
##' (see \code{\link{ncvreg}}); ignored when \code{method = "alasso"}.
##' @return An object of class \code{"polywog"}, a list containing: \describe{
##'   \item{\code{coefficients}}{the estimated coefficients.}
##'   \item{\code{lambda}}{the penalization parameter selected by
##' cross-validation.}
##'   \item{\code{fitted.values}}{the fitted mean values for each observation
##' used in fitting.}
##'   \item{\code{lmcoef}}{coefficients from an unpenalized least-squares
##' regression.}
##'   \item{\code{penwt}}{penalty weights if \code{method = "alasso"},
##' \code{NULL} otherwise.}
##'   \item{\code{formula}}{model formula.}
##'   \item{\code{degree}}{degree of polynomial basis expansion.}
##'   \item{\code{family}}{model family, \code{"gaussian"} or
##' \code{"binomial"}.}
##'   \item{\code{weights}}{observation weights if specified, \code{NULL}
##' otherwise.}
##'   \item{\code{method}}{regularization method, \code{"alasso"} or
##' \code{"scad"}.}
##'   \item{\code{penwt.method}}{estimator for penalty weights in adaptive LASSO
##' models, \code{"lm"} or \code{"glm"}; is \code{NULL} if \code{method =
##' "scad"}.}
##'   \item{\code{nfolds}}{number of cross-validation folds.}
##'   \item{\code{terms}}{the \code{\link{terms}} object used in fitting.}
##'   \item{\code{pivot}}{indices of the non-collinear columns of the full basis
##' expansion of the input variables.}
##'   \item{\code{nobs}}{number of observations used in fitting.}
##'   \item{\code{na.action}}{information on how \code{NA}s in the input data were
##' handled.}
##'   \item{\code{xlevels}}{levels of factor variables used in fitting.}
##'   \item{\code{polyTerms}}{a matrix recording how many powers of each raw
##' input term are represented in each column of the design matrix.}
##'   \item{\code{varNames}}{names of the raw input variables included in the
##' model.}
##'   \item{\code{call}}{the original function call.}
##'   \item{\code{model}}{(if requested) the model frame (see
##' \code{\link{model.frame}}).}
##'   \item{\code{X}}{(if requested) the basis-expanded design matrix.}
##'   \item{\code{y}}{(if requested) the response variable.}
##'   \item{\code{boot.matrix}}{(if \code{boot > 0}) a sparse matrix of class
##' \code{"\link[=dgCMatrix-class]{dgCMatrix}"} containing bootstrap
##' coefficients (see \code{\link{bootPolywog}}).}
##' }
##' @seealso To estimate variation via the bootstrap, see
##' \code{\link{bootPolywog}}.  To generate fitted values, see
##' \code{\link{predVals}} (and the underlying method
##' \code{\link{predict.polywog}}).  For plots, see \code{\link{plot.polywog}}.
##' The polynomial degree may be selected via cross-validation using
##' \code{\link{cv.polywog}}.
##'
##' Polynomial basis expansions of matrix inputs are computed with
##' \code{\link{polym2}}.
##'
##' Adaptive LASSO estimates are provided via \code{\link{glmnet}} and
##' \code{\link{cv.glmnet}} from the \pkg{glmnet} package.  SCAD estimates are
##' via \code{\link{cv.ncvreg}} and \code{\link{ncvreg}} in the \pkg{ncvreg}
##' package.
##' @references Brenton Kenkel and Curtis S. Signorino.  2012.  "A Method for
##' Flexible Functional Form Estimation: Bootstrapped Basis Regression with
##' Variable Selection."  Typescript, University of Rochester.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @examples
##' ## Using occupational prestige data
##' data(Prestige, package = "car")
##' Prestige <- transform(Prestige, income = income / 1000)
##'
##' ## Fit a polywog model with bootstrap iterations
##' set.seed(22)
##' fit1 <- polywog(prestige ~ education + income + type, data = Prestige,
##'                 boot = 10)
##'
##' ## Basic information
##' print(fit1)
##' summary(fit1)
##'
##' ## See how fitted values change with education holding all else fixed
##' predVals(fit1, "education", n = 10)
##'
##' ## Plot univariate relationships
##' plot(fit1)
##'
##' ## Change regularization method
##' fit2 <- update(fit1, method = "scad")
##' cbind(coef(fit1), coef(fit2))
polywog <- function(formula, data, subset, weights, na.action,
                    degree = 3,
                    family = c("gaussian", "binomial"),
                    method = c("alasso", "scad", "none"),
                    penwt.method = c("lm", "glm"),
                    boot = 0, control.boot = control.bp(),
                    .parallel = FALSE,
                    model = TRUE, X = FALSE, y = FALSE,
                    nfolds = 10, scad.maxit = 5000)
{
    cl <- match.call()
    family <- match.arg(family)
    method <- match.arg(method)
    penwt.method <- match.arg(penwt.method)
    ret.X <- X
    ret.y <- y

    ## Assemble the model frame the usual way
    formula <- as.Formula(formula)
    mf <- match(c("data", "subset", "weights", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formula
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    terms <- attr(mf, "terms")

    ## Store variable names (for use in 'print.polywog')
    varNames <- character()
    expanded <- logical()
    for (i in seq_len(length(formula)[2])) {
        tt <- terms(stats::formula(formula, rhs = i), data = mf,
                    simplify = TRUE)
        tt <- attr(tt, "term.labels")
        varNames <- c(varNames, tt)
        expanded <- c(expanded, rep(i == 1, length(tt)))
    }
    attr(varNames, "expanded") <- expanded

    ## Extract response variable and perform sanity check
    y <- model.part(formula, mf, lhs = 1, drop = TRUE)
    if (family == "binomial" && !all(y %in% c(0, 1)))
        stop("Response must be binary when family = \"binomial\"")

    ## Extract weights if any (same code as in 'lm')
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    nowt <- is.null(w)
    if (!nowt && method == "scad")
        stop("weights not allowed with method = \"scad\"")
    if (nowt)
        w <- rep(1, length(y))
    if (any(w < 0))
        stop("negative weights not allowed")

    ## Extract model matrix
    X <- makeX(formula, mf, degree)

    ## Compute linear model coefficients and eliminate singularities in X via
    ## the QR decomposition
    qx <- qr(sqrt(w) * cbind(1L, X))
    pivot <- qx$pivot[seq_len(qx$rank)]
    lmcoef <- qr.coef(qx, sqrt(w) * y)[pivot]
    qx <- NULL  # To save memory
    pivot <- pivot[-1] - 1  # Account for lack of intercept in X
    polyTerms <- attr(X, "polyTerms")[pivot, ]
    X <- X[, pivot, drop = FALSE]

    ## Compute penalty weights for adaptive lasso models
    if (method != "scad") {
        penwt <- 1 / abs(lmcoef[-1])
        if (penwt.method != "lm" && family == "binomial") {
            penwt <- penaltyWeightsBinary(X, y, w)
        } else if (penwt.method != "lm") {
            warning("For Gaussian model, forcing penwt.method = \"lm\"")
            penwt.method <- "lm"
        }
    } else {
        penwt <- NULL
        penwt.method <- NULL
    }

    ## Compute cross-validated model fit
    fitPolywog <- switch(method,
                         alasso = fitALasso,
                         scad = fitSCAD,
                         none = function(...) NULL)
    pfit <- fitPolywog(X = X, y = y, weights = w, family = family,
                       penwt = penwt, lambda = NULL, nfolds = nfolds,
                       scad.maxit = scad.maxit)

    ## Assemble object to return.  Ordering schema is as follows:
    ##   (1) Final estimates (including penalization parameter)
    ##   (2) Auxiliary quantities computed from those estimates
    ##   (3) First-stage estimates (including penalty weights)
    ##   (4) Information used in fitting (degree, observation weights, penalty
    ##   weight method, etc.)
    ##   (5) Other auxiliary information to be used for bootstrapping,
    ##   prediction, or marginal effects
    ##   (6) Other auxiliary information used for print/summary methods
    ##   (7) Optional components (model, X, y)
    ans <- list(coefficients = pfit$coef,
                lambda = pfit$lambda,
                # (2)
                fitted.values =
                if (method != "none") drop(cbind(1L, X) %*% pfit$coef) else NULL,
                # (3)
                lmcoef = lmcoef,
                penwt = penwt,
                # (4)
                formula = formula,
                degree = degree,
                family = family,
                weights = if (nowt) NULL else w,
                method = method,
                penwt.method = penwt.method,
                nfolds = nfolds,
                # (5)
                terms = terms,
                pivot = pivot,
                nobs = nrow(X),
                na.action = attr(mf, "na.action"),
                xlevels = .getXlevels(terms, mf),
                polyTerms = polyTerms,
                # (6)
                varNames = varNames,
                call = cl)
    class(ans) <- "polywog"
    if (model)
        ans$model <- mf
    if (ret.X)
        ans$X <- X
    if (ret.y)
        ans$y <- y

    ## Bootstrapping, if requested
    if (boot > 0 && method != "none") {
        ## Temporarily include model in object if necessary
        if (!model && !(ret.X && ret.y))
            ans$model <- mf

        ## Construct call to 'bootPolywog' from supplied options
        ans$boot.matrix <- do.call(bootPolywog,
                                   c(control.boot,
                                     list(model = ans, nboot = boot, .matrixOnly
                                          = TRUE, .parallel = .parallel)))

        ## Remove temporarily included model if not requested
        if (!model)
            ans$model <- NULL
    }

    return(ans)
}

##' Auxiliary for bootstrap options
##'
##' Used to pass options to \code{\link{bootPolywog}} when bootstrapping via the
##' \code{boot} option of the main \code{\link{polywog}} function.
##' @inheritParams bootPolywog
##' @return A list containing the function arguments.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
control.bp <- function(reuse.lambda = FALSE, reuse.penwt = FALSE,
                       report = FALSE, scad.maxit = 5000)
{
    list(reuse.lambda = reuse.lambda, reuse.penwt = reuse.penwt,
         report = report, scad.maxit = scad.maxit)
}

##
## Innards of the bootstrap procedure
##
bootFit <- function(X, y, weights, family, lambda, penwt, method, penwt.method,
                    nfolds, scad.maxit, pb, i)
{
    ## Calculate penalty weights unless specified
    if (method == "alasso" && is.null(penwt)) {
        if (family == "gaussian" || penwt.method == "lm") {
            penwt <- lsfit(X, y, wt = weights, intercept = TRUE)$coef
            penwt <- 1 / abs(penwt[-1])
        } else {
            penwt <- penaltyWeightsBinary(X, y, weights)
        }
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
                        reuse.penwt = FALSE, report = FALSE, scad.maxit = 5000,
                        .parallel = FALSE, .matrixOnly = FALSE)
{
    ## Load 'foreach' package for parallelization if requested
    if (.parallel && !require("foreach")) {
        .parallel <- FALSE
        warning("Must have 'foreach' package installed to use parallelization")
    }

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

    ## Bootstrap iterations
    pb <- if (report) txtProgressBar(min = 0, max = nboot) else NULL

    ## Loop over bootstrap iterations
    if (.parallel) {
        ## Loop in parallel via 'foreach'
        ans <- foreach (i = seq_len(nboot), .packages = "polywog") %dopar% {
            repeat {
                ## Ensure that the bootstrap X matrix has the same rank as the
                ## original one
                ind <- sample(seq_len(nobs), nobs, replace = TRUE)
                if (qr(cbind(1L, X[ind, , drop = FALSE]))$rank == ncf)
                    break
            }
            polywog:::bootFit(X = X[ind, , drop = FALSE], y = y[ind], weights =
                    weights[ind], family = family, lambda = lambda, penwt =
                    penwt, method = method, penwt.method = penwt.method,
                    nfolds = nfolds, scad.maxit = scad.maxit, pb = pb, i = i)
        }
    } else {
        ## Use a sequential 'for' loop
        ans <- vector("list", nboot)
        for (i in seq_len(nboot)) {
            repeat {
                ## Ensure that the bootstrap X matrix has the same rank as the
                ## original one
                ind <- sample(seq_len(nobs), nobs, replace = TRUE)
                if (qr(cbind(1L, X[ind, , drop = FALSE]))$rank == ncf)
                    break
            }
            ans[[i]] <-
                bootFit(X = X[ind, , drop = FALSE], y = y[ind], weights =
                        weights[ind], family = family, lambda = lambda, penwt =
                        penwt, method = method, penwt.method = penwt.method,
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

##' @S3method print polywog
print.polywog <- function(x, ...)
{
    ## Print the function call used to fit the model (using same code as in
    ## 'print.lm')
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    ## Extract lists of expanded and non-expanded terms
    expanded <- attr(x$varNames, "expanded")
    inTerms <- paste(x$varNames[expanded], collapse = ", ")
    cat("Variables included in polynomial expansion of degree ", x$degree, ":\n",
        sep = "")
    writeLines(strwrap(inTerms, prefix = "  "))
    if (any(!expanded)) {
        outTerms <- paste(x$varNames[!expanded], collapse = ", ")
        cat("\nVariables included linearly:\n")
        writeLines(strwrap(outTerms, prefix = "  "))
    }

    ## Print regularization method
    cat("\nRegularization method:",
        switch(x$method, alasso = "Adaptive LASSO", scad = "SCAD"))

    cat("\n\n")
    invisible(x)
}

##' @S3method vcov polywog
vcov.polywog <- function(object, ...)
{
    ncf <- length(coef(object))
    if (!is.null(object$boot.matrix)) {
        ans <- var(as.matrix(object$boot.matrix))
    } else {
        ans <- matrix(NA, nrow = ncf, ncol = ncf)
    }
    return(ans)
}

##' Summarize a fitted polywog model
##'
##' Generates a "regression table" to summarize the fitted model, including
##' coefficients along with their bootstrapped standard errors and confidence
##' intervals.  If the fitted model does not have a \code{boot.matrix} element,
##' the output will contain \code{NA}s for the standard errors, and confidence
##' intervals will not be displayed.
##' @param object a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param level width of the bootstrap confidence interval to compute for the
##' model coefficients.
##' @param prop0 logical: whether to print the proportion of bootstrap
##' iterations in which each coefficient was estimated as exactly 0.  This may
##' be informative but should \emph{not} be interpreted as a p-value.
##' @param ... other arguments, currently ignored.
##' @return An object of class \code{"summary.polywog"} whose elements are the
##' "regression table" (\code{coefficients}) and additional information from the
##' original fitted model.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method summary polywog
##' @export
##' @importMethodsFrom Matrix colMeans
summary.polywog <- function(object, level = .95, prop0 = FALSE, ...)
{
    ans <- list()

    ## Create matrix of summary results to be printed
    cf <- coef(object)
    se <- sqrt(diag(vcov(object)))
    ans$coefficients <- cbind("Estimate" = cf, "Std. Error" = se)
    if (!is.null(object$boot.matrix)) {
        q <- 0.5 - (level/2)
        interval <- apply(object$boot.matrix, 2, quantile, probs = c(q, 1-q))
        p0 <- colMeans(object$boot.matrix == 0)
        ans$coefficients <- cbind(ans$coefficients, t(interval),
                                  "Prop. 0" = if (prop0) p0)
    }

    ## Add relevant information for printing model summary
    ans$call <- object$call
    ans$degree <- object$degree
    ans$family <- object$family
    ans$method <- object$method
    ans$penwt.method <- object$penwt.method
    ans$nobs <- object$nobs
    ans$nboot <-
        if (!is.null(object$boot.matrix)) nrow(object$boot.matrix) else 0
    ans$lambda <- object$lambda

    class(ans) <- "summary.polywog"
    return(ans)
}

##' @S3method print summary.polywog
print.summary.polywog <- function(x, digits = max(3, getOption("digits") - 3),
                                  ...)
{
    ## Print the function call used to fit the model (using same code as in
    ## 'print.lm')
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    ## Print a typical 'summary'-style matrix
    cat("Coefficients:\n")
    printCoefmat(coef(x), digits = digits, ...)

    ## Print important info about the model
    cat("\nRegularization method:",
        switch(x$method, alasso = "Adaptive LASSO", scad = "SCAD"))
    if (x$method == "alasso") {
        pname <- switch(x$penwt.method,
                        lm = "inverse linear model coefficients",
                        glm = "inverse logistic regression coefficients")
        cat("\nAdaptive weights:", pname)
    }
    cat("\nNumber of observations:", x$nobs)
    cat("\nPolynomial expansion degree:", x$degree)
    cat("\nModel family:", x$family)
    cat("\nBootstrap iterations:", x$nboot)
    cat("\nPenalization parameter (lambda):", format(x$lambda, digits = digits))
    cat("\n\n")

    invisible(x)
}
