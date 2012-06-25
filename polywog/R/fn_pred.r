##' @include helpers.r
NULL

##' Predict method for polywog objects
##'
##' Generates fitted values, including bootstrap confidence intervals, for in-
##' and out-of-sample data from a fitted polywog model.
##'
##' There are some special considerations involving elements of \code{newdata}
##' that contain \code{NA}s.  If all estimated coefficients involving the
##' variable \code{z} are exactly 0 (i.e., \code{z} is completely excluded from
##' the estimated model), then it is still possible to generate a fitted value
##' for an observation where \code{z} is not recorded.  Similarly, if all
##' coefficients involving \code{z} are estimated as 0 in every bootstrap
##' iteration, it is still possible to compute confidence intervals for such an
##' observation.
##'
##' By default, \code{predict.polywog} computes fitted values and confidence
##' intervals whenever it is possible to do so.  To compute these values only
##' for fully complete cases in \code{newdata}, set \code{na.action =
##' \link{na.omit}} or \code{na.action = \link{na.exclude}}.  \code{na.omit}
##' will return an object whose dimension equals the number of complete cases;
##' \code{na.exclude} gives the same dimension as \code{newdata}, padded with
##' \code{NA}s where appropriate.
##' @param object a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param newdata an optional data frame containing observations for which
##' fitted values should be computed.  If not specified, fitted values are
##' generated for the data used to fit the model.
##' @param type specifies whether the fitted values should be generated on the
##' link scale (\eqn{X \beta}) or in terms of the expected value of the response
##' variable.  These only differ for binomial family models.
##' @param interval logical: whether to calculate bootstrap confidence intervals
##' for each fitted value.
##' @param level confidence level for the intervals.
##' @param bag logical: whether to use "bootstrap aggregation" to generate the
##' main fitted values (if \code{FALSE}, they are calculated from the main model
##' fit).
##' @param na.action a function specifying what to do with observations in
##' \code{newdata} containing \code{NA}s (default \code{\link{na.pass}}).  See
##' "Details".
##' @param ... additional arguments, currently ignored.
##' @return If \code{interval = TRUE}, a matrix containing each fitted value and
##' its confidence interval.  Otherwise, a vector containing the fitted values.
##' @seealso For more user-friendly generation of fitted values, see
##' \code{\link{predVals}}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @method predict polywog
##' @export
predict.polywog <- function(object, newdata,
                            type = c("link", "response"),
                            interval = FALSE, level = .95,
                            bag = FALSE,
                            na.action = na.pass, ...)
{
    type <- match.arg(type)
    transform <- (type == "response" && object$family == "binomial")

    ## Check for nothing that depends on bootstrap results if 'object' does not
    ## have a 'boot.matrix' element
    if (is.null(object$boot.matrix) && (interval || bag))
    {
        interval <- bag <- FALSE
        warning("Options 'interval' and 'bag' not available for models without a 'boot.matrix' element")
    }

    ## Setup largely adapted from predict.lm() code; the bits relating to
    ## 'newdata' being a model frame are adapted from mgcv::predict.gam()
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        ## Use original model matrix if available
        X <- object$X
        if (is.null(X)) {
            if (is.null(object$model))
                stop("Fitted object must contain either 'model' or 'X' to use predict.polywog without specifying 'newdata'; re-run polywog with \"model = TRUE\"")
            X <- makeX(object$formula, object$model, object$degree)
            X <- X[, object$pivot, drop = FALSE]
        }
        nd.is.mf <- FALSE
    } else if (is.data.frame(newdata) && !is.null(attr(newdata, "terms"))) {
        ## 'newdata' is a model frame -- this case must be treated separately,
        ## or else predVals() and margEff.polywog() won't work when the original
        ## model formula contains transformations of the original inputs

        X <- makeX(object$formula, newdata, object$degree, na.ok = TRUE)
        X <- X[, object$pivot, drop = FALSE]
        nd.is.mf <- TRUE
    } else {
        ## Construct model frame from 'newdata'
        Terms <- delete.response(tt)
        mf <- model.frame(Terms, newdata, na.action = na.action, xlev =
                          object$xlevels)

        ## Check validity of 'newdata' (covariate types same as in fitted model)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, mf)

        ## Construct polynomial-expanded model matrix and remove collinear terms
        X <- makeX(object$formula, mf, object$degree, na.ok = TRUE)
        X <- X[, object$pivot, drop = FALSE]
        nd.is.mf <- FALSE
    }
    X <- cbind("(Intercept)" = 1L, X)

    ## Generate predicted values on link scale
    ##
    ## Special things are done with "zeroed out" coefficients to ensure that
    ## predictions can be returned for variables that are NA if they are not
    ## "included" in the model
    inMainModel <- coef(object) != 0
    if (bag || interval) {
        allPreds <- apply(object$boot.matrix, 1, function(b) {
            inBootModel <- b != 0
            X[, inBootModel, drop = FALSE] %*% b[inBootModel]
        })
        if (transform)
            allPreds <- plogis(allPreds)
    }
    if (bag) {
        pred <- rowMeans(allPreds)
    } else if (transform) {
        pred <- plogis(drop(X[, inMainModel, drop = FALSE] %*%
                            coef(object)[inMainModel]))
    } else {
        pred <- drop(X[, inMainModel, drop = FALSE] %*%
                     coef(object)[inMainModel])
    }

    ## Confidence intervals
    if (interval) {
        if (level < 0 || level > 1)
            stop("'level' must be between 0 and 1")
        q <- 0.5 - (level/2)

        ## The 'tryCatch' statements are necessary to ensure the whole operation
        ## doesn't fail due to an NA in a single row
        lwr <- apply(allPreds, 1,
                     function(x) tryCatch(quantile(x, probs = q),
                                          error = function(e) NA))
        upr <- apply(allPreds, 1,
                     function(x) tryCatch(quantile(x, probs = 1-q),
                                          error = function(e) NA))
        pred <- cbind(fit = pred, lwr = lwr, upr = upr)
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
    
    return(pred)
}

##
## Match supplied character string to variable names in data frame, and warn
## about no match or multiple partial matches.  Used in 'predVals', not meant to
## be called directly by users.
##
## ARGUMENTS:
##   vars: character string of supplied names
##   names: character string of names in data
##
## RETURN:
##   Indices of elements in 'names' matching elements in 'vars
##
getXcols <- function(vars, names)
{
    ## Figure out the index of 'names' corresponding to each value in 'vars' (NA
    ## for no match, 0 for multiple partial matches)
    vcols <- charmatch(vars, names)
    if (any(is.na(vcols))) {
        stop("The following variable names were not found in the fitted model: ",
             paste(vars[is.na(vcols)], collapse = ", "))
    } else if (any(vcols == 0)) {
        stop("The following names partially match mutliple variables in the fitted model: ",
             paste(vars[vcols==0], collapse = ", "))
    }

    return(vcols)
}

##
## Create the varying part of a fitted value profile.  Used in 'predVals', not
## meant to be called directly by users.
##
## ARGUMENTS:
##   data: data frame containing raw data
##   cols: columns of varying variables
##   xlims: named list specifying limits to use for continuous variables
##     (default is to use observed range)
##   n: number of grid points for continuous variables
##
## RETURN:
##   Data frame containing grid values for specified variables
##
getXvals <- function(data, cols, xlims, n)
{
    ans <- list()
    for (i in seq_along(cols)) {
        xcol <- cols[i]
        xname <- names(data)[xcol]
        x <- data[, xcol]

        ## Construct the sequence of values of the variable to evaluate at
        if (all(unique(x) %in% c(0, 1))) {
            ## Binary variable: use 0 and 1
            xs <- c(0, 1)
        } else if (is.numeric(x)) {
            ## Continuous variable: make a grid of size 'n' within the specified
            ## limits (or observed range if none specified)
            xlim <- xlims[[xname]]
            if (is.null(xlim))
                xlim <- range(x)
            xs <- seq(xlim[1], xlim[2], length.out = n)
        } else if (is.factor(x)) {
            ## Factor variable: use all levels
            xs <- factor(levels(x), levels = levels(x))
        } else if (is.logical(x)) {
            ## Logical variable: use TRUE and FALSE
            xs <- c(FALSE, TRUE)
        }

        ans[[i]] <- xs
        names(ans)[i] <- xname
    }

    ## Create all possible combinations of the specified variables
    ans <- expand.grid(ans)
    return(ans)
}

## makeProfile: same as in the 'games' package (see games/R/predProbs.r)
##
## INPUT:
## x: a data frame
## ...: expressions (see below)
##
## RETURN:
## a one-row data frame
##
## The function takes a data frame 'x' and returns a one-row data frame with the
## same variables, containing a "typical observation" profile from x.  The
## default is to take the mean of numeric variables, the median of ordered
## and binary variables, and the mode of categorical variables.
##
## These defaults can be overridden by passing expressions to "...".  For
## example, to set a variable 'z' to its 25th percentile, use:
##     makeProfile(x, z = quantile(z, 0.25))
## To set 'z' to equal 0.3 and 'w' to equal its maximum, use:
##     makeProfile(x, z = 0.3, w = max(w))
##
## This function is used in 'predProbs' and is not meant to be called by users.
## It is analogous to the 'setx' function in the Zelig package.
##
##' @importFrom games Mode
makeProfile <- function(x, ...)
{
    cl <- match.call(expand.dots = FALSE)

    ## get the first row of the data frame (and continue to store as data frame
    ## to allow for different data types, preserve factor levels, etc)
    ans <- x[1, ]

    ## loop over each variable in the data frame
    for (i in seq_len(ncol(x))) {
        xvar <- x[, i]

        ## use mean for numeric variables, median for ordered or dummy
        ## variables, mode for categorical variables
        isDummy <- all(unique(xvar) %in% c(0, 1))
        if (is.numeric(xvar) && !isDummy) {
            ans[, i] <- mean(xvar)
        } else if (is.ordered(xvar)) {
            ans[, i] <- median(xvar)
        } else if (isDummy) {
            ans[, i] <- median(xvar)
        } else {
            ans[, i] <- Mode(xvar)
        }
    }

    if ("..." %in% names(cl)) {        
        ## takes the expressions fed to "...", evaluates them within the
        ## supplied data frame, and returns them as a vector.
        ##
        ## e.g., if the call is makeProfile(x, foo = median(foo), bar =
        ## quantile(bar, .25)), this will return a vector with named elements
        ## "foo" and "bar".
        ##
        ## the eval(substitute()) business is to ensure that these aren't
        ## evaluated in the global environment rather than the data frame "x",
        ## which would most likely throw an error when "foo" and "bar" weren't
        ## found in the workspace, or possibly obtain the wrong values for
        ## these.
        dots <- eval(substitute(list(...)), x)
        dots <- unlist(lapply(dots, unname))

        toReplace <- match(names(dots), names(x), 0L)
        ans[toReplace] <- dots[toReplace != 0L]
    }

    ## ensures that choices in "..." expressed as characters (e.g., foo = "a"
    ## for a factor foo with levels a, b, c) don't wind up turning factor
    ## variables into characters
    fvars <- sapply(x, inherits, what = "factor")
    for (i in seq_len(ncol(x))) {
        if (fvars[i]) {
            ans[, i] <- factor(ans[, i], levels = levels(x[, i]))
        }
    }
    
    return(ans)
}


##
## Predicted values for "typical value" profiles from polywog models
##
## ARGUMENTS:
##   model: a fitted model of class "polywog"
##   boot.matrix: bootstrap coefficients
##   xvars: character string naming the covariates to be varied in the profile
##     (partial matches ok)
##   xlims: optional named list (partial matches not allowed) specifying the
##     limits to use for each continuous variable; the observed range will be
##     used for those not specified
##   n: number of grid points to use for continuous variables
##   interval, level, bag: same as in 'predict.polywog'
##   maxrows: upper limit on number of rows in the profile constructed (this is
##     to prevent from accidentally trying to calculate 1 million fitted
##     values after specifying three continuous variables)
##
## RETURN:
##   Data frame with fitted values ('fit'), confidence interval if specified
##   ('lwr' and 'upr'), and the profile data.
##

##' Easy computation of fitted values
##'
##' User-friendly generation of fitted values and their confidence intervals
##' from models of class \code{"polywog"}, using \code{\link{predict.polywog}}
##' as a backend.
##'
##' \code{predVals} allows users to examine the estimated effects of input
##' variables on the expected outcome using the coefficients returned by
##' \code{\link{polywog}}.  The procedure is designed so that, for a preliminary
##' analysis, the user can simply specify the fitted model and the independent
##' variable of interest, and quickly obtain predicted values.  However, it is
##' flexible enough to allow for finely tuned analysis as well.  The function is
##' very similar to \code{\link[games]{predProbs}} in the \pkg{games} package.
##' 
##' The procedure works by varying \code{xvars}, the variables of interest,
##' across their observed ranges (or those specified by the user in
##' \code{xlims}) while holding all other independent variables in the model
##' fixed.  The profile created by default is as follows (the same defaults as
##' in the \code{\link[Zelig]{sim}} function in the \pkg{Zelig} package):
##' \itemize{
##'   \item numeric, non-binary variables are fixed at their means
##'   \item \code{\link{ordered}} and binary variables are fixed at their
##' medians
##'   \item all others are fixed at their modes (see \code{\link[games]{Mode}})
##' }
##' However, it is possible to override these defaults for any or all
##' variables.  For example, to set a variable named \code{polity} to its lower
##' quartile, call \code{predVals} with the argument \code{polity =
##' quantile(polity, 0.25)}.  To set a factor variable to a particular level,
##' provide the name of the level as a character string (in quotes).  See the
##' examples below for illustrations of this functionality.
##'
##' All confidence intervals are generated via the bootstrap.  If \code{model}
##' does not have a \code{boot.matrix} element (see \code{\link{bootPolywog}}),
##' confidence intervals will not be computed.
##' @param model a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param xvars a character vector containing names of raw input variables
##' (from \code{model$varNames}).  Partial matches are allowed.
##' @param xlims named list of limits for the evaluation grid for each
##' continuous variable in \code{xvars}.  If not given, the variable's observed
##' range is used.
##' @param n number of grid points at which to evaluate each continuous variable
##' in \code{xvars}.
##' @param interval logical: whether to compute bootstrap confidence intervals
##' for each fitted value.
##' @param level confidence level for the intervals.
##' @param bag logical: whether to use "bootstrap aggregation" to generate the
##' main fitted values (if \code{FALSE}, they are calculated from the main model
##' fit).
##' @param maxrows maximum number of rows of output.  Used to prevent accidental
##' memory overruns when \code{xvars} contains more than two continuous
##' variables.
##' @param ... used to set values for the variables other than \code{xvars} in
##' the profile of observations.  See "Details" below.
##' @return A data frame containing the fitted values and confidence intervals
##' (if requested) for each combination of covariate values.
##' @seealso \code{\link{predict.polywog}} for more flexible (but less
##' user-friendly) computation of fitted values.  \code{\link{plot.polywog}} for
##' plotting fitted values and their confidence intervals.
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
##' ## Predicted prestige across occupational categories
##' predVals(fit1, "type")
##' predVals(fit1, "type", income = quantile(income, 0.25))
##'
##' ## Predicted prestige by education
##' predVals(fit1, "education", n = 10)
##' predVals(fit1, "education", n = 10, income = quantile(income, 0.25))
predVals <- function(model, xvars, xlims = list(), n = 100, interval = TRUE,
                     level = .95, bag = FALSE, maxrows = 10000, ...)
{
    ## Check for nothing that depends on bootstrap results if 'boot.matrix' is
    ## not supplied
    if (is.null(model$boot.matrix) && (interval || bag))
    {
        interval <- bag <- FALSE
        warning("Options 'interval' and 'bag' not available for models without a 'boot.matrix' element")
    }

    ## Extract original model data and remove dependent variable
    dat <- model$model
    if (is.null(dat))
        stop("Fitted object must contain 'model' to use predVals; re-run polywog with \"model = TRUE\"")
    dat <- model.frame(delete.response(terms(model)), data = dat, na.action =
                       na.omit, xlev = model$xlevels)

    ## Create the varying part of the fitted value profile
    xc <- getXcols(xvars, names(dat))
    xv <- getXvals(dat, xc, xlims, n)
    if (nrow(xv) > maxrows) {
        stop("Profile generated is too large; re-run with lower n or maxrows >= ",
             nrow(xv), " to continue")
    }

    ## Create the fixed part of the profile and combine with varying part
    xf <- makeProfile(dat, ...)
    x <- xf[rep(1, nrow(xv)), , drop = FALSE]
    rownames(x) <- seq_len(nrow(x))
    x[, xc] <- xv

    ## Calculate fitted values and combine with data frame
    fit <- predict(model, newdata = x, type = "response", interval = interval,
                   level = level, bag = bag)
    ans <- structure(cbind(fit, x),
                     interval = interval,
                     xvars = xvars,
                     xcol = xc + ifelse(length(dim(fit)), 3, 1))
    return(ans)
}
