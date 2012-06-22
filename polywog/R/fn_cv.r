##' @include helpers.r
##' @include fitters.r
##' @include polywog.r
NULL

##' Cross-validation for polywog models
##'
##' Uses k-fold cross-validation to select the polynomial degree and
##' penalization parameter for a \code{\link{polywog}} model.
##' @param formula model formula specifying the response and input variables.
##' @param ... other arguments to be passed to \code{\link{polywog}}.  Arguments
##' controlling the bootstrap will be ignored.
##' @param method variable selection method: \code{"alasso"} (default) for
##' adaptive LASSO or \code{"scad"} for SCAD.
##' @param model logical: whether to include the model frame in the returned
##' object (default \code{TRUE}).  It may be problematic to run
##' \code{\link{predVals}} or \code{\link{bootPolywog}} on a polywog object fit
##' with \code{model = FALSE}.
##' @param X logical: whether to include the polynomial-expanded design matrix
##' in the returned object (default \code{FALSE}).
##' @param y logical: whether to include the response variable in the returned
##' object (default \code{FALSE}).
##' @param degrees.cv vector of polynomial degrees to examine via
##' cross-validation.
##' @param nfolds number of folds to use in cross-validation to select the
##' penalization parameter.
##' @param scad.maxit maximum number of iterations when \code{method = "scad"}
##' (see \code{\link{ncvreg}}); ignored when \code{method = "alasso"}.
##' @return An object of class \code{"cv.polywog"}, a list containing:
##' \describe{
##'   \item{\code{results}}{A table of each degree tested, the optimal
##' penalization parameter \eqn{\lambda} for that degree, and its
##' cross-validation error.}
##'   \item{\code{degree.min}}{The polynomial degree giving the lowest
##' cross-validation error.}
##'   \item{\code{polywog.fit}}{A \code{\link{polywog}} model, fit at the
##' polynomial degree giving the lowest cross-validation error.}
##' }
##'
##' Because the returned object contains the fitted polywog model for the
##' optimal degree, no additional runs of \code{\link{polywog}} are necessary to
##' estimate coefficients or the penalization parameter \eqn{\lambda}.  However,
##' bootstrap coefficients must be obtained by running \code{\link{bootPolywog}}
##' on the \code{"polywog.fit"} element of the returned object, as in the
##' examples below.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @examples
##' ## Using occupational prestige data
##' data(Prestige, package = "car")
##' Prestige <- transform(Prestige, income = income / 1000)
##'
##' ## Examine degrees 1 through 4
##' set.seed(39)
##' cv1 <- cv.polywog(prestige ~ education + income + type, data = Prestige,
##'                   degrees.cv = 1:4, nfolds = 10)
##'
##' print(cv1)
cv.polywog <- function(formula,
                       ...,
                       method = c("alasso", "scad"),
                       model = TRUE, X = FALSE, y = FALSE,
                       degrees.cv = 1:3,
                       nfolds = 10, scad.maxit = 5000)
{
    cl <- match.call()
    method <- match.arg(method)

    cvfit <- switch(method, alasso = cvALasso, scad = cvSCAD)

    ## Loop over each degree to be considered
    ans <- calls <- vector("list", length(degrees.cv))
    for (i in seq_along(degrees.cv)) {
        d <- degrees.cv[i]

        ## Use the 'polywog' function with method = "none" to extract the model
        ## matrix, penalty weights, and any other relevant info.  Need to call
        ## the function from here or else the arguments in '...' get messed up
        fit <- match(names(formals(polywog)), names(cl), 0L)
        fit <- cl[c(1L, fit)]
        fit$degree <- d
        fit$model <- FALSE
        fit$X <- TRUE
        fit$y <- TRUE
        fit$method <- "none"
        fit[[1]] <- as.name("polywog")
        calls[[i]] <- fit  # save so can call the best one again later
        fit <- eval(fit, parent.frame())

        ## In the first iteration, randomly assign cross-validation folds
        ## (this is done inside the loop so we can access the number of
        ## observations from the fitted polywog object)
        if (i == 1L) {
            n <- fit$nobs
            foldid <- sample(rep(seq_len(nfolds), length.out = n))
        }

        ## Cross-validate and store results
        ans[[i]] <- cvfit(X = fit$X, y = fit$y,
                          weights = if (!is.null(fit$weights)) fit$weights
                          else rep(1, n),
                          family = fit$family, penwt = fit$penwt,
                          nfolds = nfolds, foldid = foldid,
                          scad.maxit = scad.maxit)

        ## Save memory
        rm(fit)
    }

    ## Table of minimal cross-validation error by degrees
    tab <- t(sapply(ans, "[", c("lambda", "cverr")))
    tab <- cbind(degrees.cv, tab)
    colnames(tab) <- c("degree", "lambda.min", "cverr")
    best <- which.min(tab[, "cverr"])

    ## Create the fitted polywog object from the degree with minimal
    ## cross-validation error.  Could have saved each one created within the
    ## loop and then chosen, but that would add a lot of memory use for little
    ## computational speedup.
    polywog.fit <- calls[[best]]
    polywog.fit$model <- model
    polywog.fit$X <- X
    polywog.fit$y <- y
    polywog.fit$model <- model
    polywog.fit$degree <- degrees.cv[best]
    polywog.fit <- eval(polywog.fit, parent.frame())
    polywog.fit$coefficients <- ans[[best]]$coef
    polywog.fit$lambda <- ans[[best]]$lambda
    polywog.fit$method <- method
    polywog.fit$call$method <- method
    if (method == "scad") {
        polywog.fit$penwt <- NULL
        polywog.fit$penwt.method <- NULL
    }

    ans <- list(results = tab, degree.min = best, polywog.fit = polywog.fit)
    class(ans) <- "cv.polywog"
    ans
}

