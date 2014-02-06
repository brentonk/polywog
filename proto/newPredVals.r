##############################################################################
### newPredVals.r
###
### Modified version of the predVals() function in the polywog package, which
### uses the "observed value" approach of Hanmer and Kalkan (AJPS 2012).
### This is a stopgap until the new version can be incorporated into the
### package.
###
### Brenton Kenkel and Curt Signorino
##############################################################################

library(foreach)
library(Rcpp)

sourceCpp("newPredVals.cpp")

newPredict <- function(object, newdata,
                       type = c("link", "response"),
                       interval = FALSE, level = .95,
                       bag = FALSE,
                       na.action = na.pass, ...)
{
    ## TODO:
    ##   o Implement bootstrap aggregation
    ##   o Implement response transformation for logit
    ##   o See whether the "nd.is.mf" check is necessary

    type <- match.arg(type)
    transform <- (type == "response" && object$family == "binomial")

    ## Check for nothing that depends on bootstrap results if 'object' does
    ## not have a 'boot.matrix' element
    if (is.null(object$boot.matrix) && (interval || bag))
    {
        interval <- bag <- FALSE
        warning("Options 'interval' and 'bag' not available for models without a 'boot.matrix' element")
    }

    ## Setup largely adapted from predict.lm() code; the bits relating to
    ## 'newdata' being a model frame are adapted from mgcv::predict.gam()
    X.exists <- FALSE
    if (missing(newdata) || is.null(newdata)) {
        ## Use original model matrix if available
        X <- object$X
        if (is.null(X)) {
            if (is.null(object$model))
                stop("Fitted object must contain either 'model' or 'X' to use predict.polywog without specifying 'newdata'; re-run polywog with \"model = TRUE\"")
            mf <- object$model
        } else {
            X.exists <- TRUE
        }
        nd.is.mf <- FALSE
    } else if (is.data.frame(newdata) && !is.null(attr(newdata, "terms"))) {
        ## 'newdata' is a model frame -- this case must be treated separately,
        ## or else predVals() and margEff.polywog() won't work when the
        ## original model formula contains transformations of the original
        ## inputs

        mf <- newdata
        nd.is.mf <- TRUE
    } else {
        ## Construct model frame from 'newdata'
        Terms <- delete.response(terms(object))
        mf <- model.frame(Terms, newdata, na.action = na.action, xlev =
                          object$xlevels)

        ## Check validity of 'newdata' (covariate types same as in fitted
        ## model)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, mf)

        nd.is.mf <- FALSE
    }

    ## Compute the model matrix
    if (!X.exists) {
        ff <- polywog:::removeIntercepts(object$formula, mf)
        X <- model.matrix(ff, data = mf, rhs = 1)
        if (length(object$formula)[2] > 1)
            X <- cbind(X, model.matrix(ff, data = mf, rhs = 2))
    }

    pred <- newPredictPolywogC(X = X, poly_terms = object$polyTerms,
                               coef = list(main = coef(object),
                               boot = if (interval) t(object$boot.matrix)),
                               avg = FALSE, interval = interval,
                               level = level)
    if (interval) {
        pred <- do.call(cbind, pred)
    } else {
        pred <- pred$pred
    }

    ## If just computing in-sample fits, ensure conformity with the original
    ## 'na.action' (i.e. padding when na.action = na.exclude, as in
    ## 'predict.lm')
    if (missing(newdata) || is.null(newdata)) {
        pred <- napredict(object$na.action, pred)
    } else if (nd.is.mf) {
        pred <- napredict(attr(newdata, "na.action"), pred)
    } else {
        pred <- napredict(attr(mf, "na.action"), pred)
    }

    pred
}

predValsR <- function(model, data = model$model,
                      xvars, xlims = list(), n = 50,
                      interval = TRUE, level = .95,
                      report = FALSE, .parallel = FALSE)
{
    ## Can't form a confidence interval if no bootstrap results
    if (is.null(model$boot.matrix) && interval) {
        interval <- FALSE
        warning("Option 'interval' not available for models without a 'boot.matrix' element")
    }

    ## Can't print a status bar in parallel
    if (.parallel)
        report <- FALSE

    ## Calculate grid of covariate values to examine
    xc <- polywog:::getXcols(xvars, names(data))
    xv <- polywog:::getXvals(data, xc, xlims, n)

    ## Loop through the grid of covariate values
    `%dofn%` <- if (.parallel) `%dopar%` else `%do%`
    if (report)
        pb <- txtProgressBar(min = 0, max = nrow(xv))
    ans <- foreach (i = seq_len(nrow(xv))) %dofn% {
        ## Replace the actual values of the selected covariates in the data
        ## with the i'th row of the grid, while keeping everything else at its
        ## true values
        data[, names(xv)] <- xv[i, ]

        ## Compute the average predicted response within the dataset when the
        ## given covariates are set at these values
        avg <- mean(predict(model, newdata = data, type = "response"))

        ## Compute confidence interval from bootstrap coefficients (this needs
        ## to be done "by hand" since we can't extract the individual
        ## bootstrap predictions via predict.polywog())
        if (interval) {
            ## Compute polynomial basis expansion of 'data'
            X <- polywog:::makeX(model$formula, data, model$degree,
                                 na.ok = TRUE)
            X <- X[, model$pivot, drop = FALSE]
            X <- cbind("(Intercept)" = 1L, X)

            ## Compute N[data] x B matrix of bootstrap predicted values
            pred <- X %*% t(model$boot.matrix)
            if (model$family == "binomial")
                pred <- plogis(pred)

            ## Compute the B bootstrap averages and find the relevant order
            ## statistics
            pred <- colMeans(pred)  # each column is one bootstrap iteration
            lwr <- unname(quantile(pred, probs = 0.5 - level/2))
            upr <- unname(quantile(pred, probs = 0.5 + level/2))
        } else {
            lwr <- upr <- NA
        }

        if (report)
            setTxtProgressBar(pb, i)

        c(avg = avg, lwr = lwr, upr = upr)
    }
    ans <- do.call(rbind, ans)
    ans <- data.frame(cbind(xv, ans))

    if (!interval)
        ans$lwr <- ans$upr <- NULL

    ans
}

