##' @include helpers.r
##' @include fn_pred.r
NULL

##' @method margEff polywog
##' @export
margEff.polywog <- function(object, xvar = NULL, drop = FALSE, ...)
{
    ## If no variable specified, apply to each variable in the model
    if (is.null(xvar)) {
        ans <- lapply(object$varNames,
                      function(x) margEff(object, xvar = x, drop = drop))
        names(ans) <- object$varNames
        class(ans) <- "margEff.polywog"
        return(ans)
    }

    ## Extract original model data
    mf <- object$model
    if (is.null(mf))
        stop("Fitted object must contain 'model' to use margEff; re-run polywog with \"model = TRUE\"")

    ## Extract the variable of interest and determine its type
    xc <- getXcols(xvar, names(mf))
    xvar <- names(mf)[xc]  # replace with completed name
    x <- mf[, xc]
    xbinary <- is.numeric(x) && all(x %in% c(0, 1))

    if (is.numeric(x) && !xbinary) {
        ##-- CONTINUOUS VARIABLE --##

        ## Extract the raw/un-transformed model matrix
        formula <- removeIntercepts(object$formula, mf)
        Xraw <- model.matrix(formula, data = mf, rhs = 1)
        if (length(object$formula)[2] > 1)
            Xraw <- cbind(Xraw, model.matrix(formula, data = mf, rhs = 2))

        xc <- getXcols(xvar, colnames(Xraw))
        ncoef <- length(coef(object)) - 1  # Not considering intercept
        pt <- object$polyTerms
        ans <- rep(0, nrow(Xraw))

        ## Loop over each coefficient
        for (i in seq_len(ncoef)) {
            ## Move on if the term doesn't include the variable of interest or
            ## if the term's coefficient is zero
            if (pt[i, xc] == 0 || coef(object)[i+1] == 0)
                next

            ## Compute appropriate powers of input variables
            powers <- pt[i, ]
            powers[xc] <- powers[xc] - 1
            Xpow <- sweep(Xraw, 2, powers, "^")

            ## Compute marginal effect with respect to this coefficient
            ans <- ans + (powers[xc] + 1) * coef(object)[i+1] * .rowProds(Xpow)
        }
    } else {
        ##-- DISCRETE VARIABLE --##

        ## Switch to determine the "levels"
        if (is.factor(x)) {
            levs <- levels(x)
        } else if (is.character(x)) {
            warning("character variables should be converted to factors before fitting with polywog(); errors or unpredictable behavior may result")
            levs <- levels(as.factor(x))
        } else if (is.logical(x)) {
            levs <- c(TRUE, FALSE)
        } else if (xbinary) {
            levs <- c(1, 0)
        }

        ## Run predict() on the original data, setting the selected variable to
        ## each of its potential levels
        ans <- vector("list", length(levs))
        for (i in seq_along(levs)) {
            newmf <- mf
            newmf[, xc][] <- levs[i]

            ## Need to ensure that the "na.action" attribute of 'newdata' is
            ## NULL, or else padding may result if the model was originally run
            ## with na.exclude
            ans[[i]] <- predict(object,
                                newdata = structure(newmf, na.action = NULL),
                                type = "response")
        }
        ans <- do.call(cbind, ans)
        colnames(ans) <- levs

        ## Take first differences with respect to the reference category (last
        ## level)
        ans <- ans - ans[, ncol(ans)]
        ans <- ans[, -ncol(ans), drop = FALSE]  # Remove column of zeroes
        if (drop)
            ans <- base::drop(ans)
        attr(ans, "levels") <- levs
    }

    attr(ans, "xvar") <- xvar
    class(ans) <- "margEff.polywog"
    ans
}

##' @S3method print margEff.polywog
print.margEff.polywog <- function(x, ...)
{
    ## Print minimal information
    if (is.list(x)) {
        cat("Marginal effects of all variables in a fitted model of class \"polywog\"\n")
    } else {
        cat("Marginal effect of", attr(x, "xvar"), "in a fitted model of class \"polywog\"\n")
    }
    invisible(x)
}


## Convenience function to make a data frame out of a list of marginal effects
MEtoDF <- function(object)
{
    ## Figure out which entries are first differences
    whichFactors <- sapply(object, function(x) !is.null(attr(x, "levels")))

    ## Make intelligible names for first difference entries
    for (i in seq_len(sum(whichFactors))) {
        j <- which(whichFactors)[i]
        object[[j]] <- as.matrix(object[[j]])
        levs <- attr(object[[j]], "levels")
        xvar <- attr(object[[j]], "xvar")
        colnames(object[[j]]) <-
            paste(xvar, "=", levs[-length(levs)], " (vs. ",
                  levs[length(levs)], ")", sep = "")
    }

    as.data.frame(do.call(cbind, object))
}

##' @S3method summary margEff.polywog
summary.margEff.polywog <- function(object, probs = seq(0, 1, by = 0.25), ...)
{
    if (!is.list(object)) {
        object <- list(object)
        names(object) <- attr(object[[1]], "xvar")
    }

    object <- MEtoDF(object)
    object <- sapply(object,
                     function(x) c("Mean" = mean(x), "SD" = sd(x),
                                   quantile(x, probs = probs)))
    object <- t(object)
    class(object) <- "summary.margEff.polywog"
    object
}

##' @S3method print summary.margEff.polywog
print.summary.margEff.polywog <- function(x,
                                          digits = max(3, getOption("digits") - 3),
                                          ...)
{
    printCoefmat(zapsmall(x), digits = digits, ...)
    invisible(x)
}

##' @S3method plot margEff.polywog
plot.margEff.polywog <- function(x, ...)
{
    if (!is.list(x)) {
        x <- list(x)
        names(x) <- attr(x[[1]], "xvar")
    }

    x <- MEtoDF(x)

    ## Obtain dimensions
    mfrow <- ceiling(sqrt(length(x)))
    mfcol <- ceiling(length(x) / mfrow)
    op <- par(mfrow = c(mfrow, mfcol))
    on.exit(par(op))

    ## Plot each density plot
    for (i in seq_along(x))
        plot(density(x[[i]]), main = names(x)[i], ...)

    invisible(x)
}

## Bug-fixed version of matrixStats::rowProds (which was released under an
## Artistic-2.0 license, which is compatible with the GPL and hence can be
## redistributed with modification here)

##' @importFrom matrixStats rowSums
.rowProds <- function (x, ...) 
{
    s <- (x == 0)
    s <- rowSums(s)
    ok <- (s == 0)
    rm(s)
    y <- vector(mode(x), nrow(x))
    x <- x[ok, , drop = FALSE]
    s <- (x < 0)
    s <- rowSums(s)
    s <- (s%%2)
    s <- c(+1, -1)[s + 1]
    x <- abs(x)
    x <- log(x)
    x <- rowSums(x, ...)
    x <- exp(x)
    x <- s * x
    y[ok] <- x
    rm(ok, s, x)
    y
}
