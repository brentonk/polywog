##############################################################################
### newPredVals.r
###
### Modified version of the predVals() function in the polywog package, which
### uses the "observed value" approach of Hanmer and Kalkan (AJPS 2012).
### This is a stopgap until the new version can be incorporated into the
### package.
###
### Brenton Kenkel and Curt Signorino
### last updated 2013-01-14 (bjk)
##############################################################################

library(foreach)

newPredVals <- function(model, data = model$model,
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

