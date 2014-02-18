preplotFromPick <- function(x, pick, ...)
{
    pp <- predVals(x, xvars = pick, ...)
    ans <- structure(pp, class = c("preplot.polywog", "data.frame"))
    return(ans)
}

##' Univariate and bivariate fitted value plots
##'
##' Generates plots of the relationship between input variables and the expected
##' value of the outcome, using \code{\link{predVals}} as a backend.
##'
##' By default, a univariate plot generated by \code{plot.polywog} shows the
##' relationship between the selected input variable and the expected outcome
##' while holding all other covariates at "central" values (as in
##' \code{\link{predVals}}).  The values that the other variables are held out
##' can be changed by supplying additional arguments to \code{...}, as in the
##' examples below.
##'
##' Similarly, a bivariate plot shows the relationship between two input
##' variables and the expected outcome while holding all else fixed.  If either
##' variable is binary or categorical, the plot will show the relationship
##' between one variable and the expected outcome across each value/level of the
##' other.
##' @param x a fitted model of class \code{"polywog"}, typically the output
##' of \code{\link{polywog}}.
##' @param which selection of variables to plot: a character vector containing
##' one or two names of raw input variables (see \code{x$varNames}).  May also
##' be a numeric vector corresponding to indices of \code{x$varNames}.
##' If \code{which = NULL}, a plot of each individual term will be generated.
##' @param ask logical: whether to display an interactive menu of terms to
##' select.
##' @param auto.set.par logical: whether to temporarily change the graphics
##' parameters so that multiple plots are displayed in one window (e.g., each
##' univariate plot when \code{which = NULL}).
##' @param interval logical: whether to display bootstrap confidence intervals
##' around each fitted value.  Not available for bivariate plots unless
##' \code{FUN3d = "persp3d"}.
##' @param level confidence level for the intervals.
##' @param bag logical: whether to use "bootstrap aggregation" to generate the
##' main fitted values (if \code{FALSE}, they are calculated from the main model
##' fit).
##' @param FUN3D which plotting function to use to generate bivariate plots.
##' Valid options include \code{"\link{contour}"} (the default) and
##' \code{"\link{filled.contour}"}; \code{"\link[lattice]{wireframe}"}, which
##' requires the \pkg{lattice} package; and \code{"\link[rgl]{persp3d}"}, which
##' requires the \pkg{rgl} package.
##' @param control.plot list of arguments to be passed to the underlying
##' plotting functions (e.g., axis labels and limits).
##' @param ... additional arguments to be passed to \code{\link{predVals}}.
##' @return An object of class \code{preplot.polywog}, invisibly.  This is a
##' data frame generated by \code{\link{predVals}} that contains all information
##' used in plotting.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @importFrom stringr str_split
##' @method plot polywog
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
##' ## All univariate relationships
##' plot(fit1)
##'
##' ## Predicted prestige across occupational categories
##' plot(fit1, which = "type",
##' control.plot = list(xlab = "occupational category"))
##'
##' ## Predicted prestige by education across occupational categories
##' plot(fit1, which = c("education", "type"))
##'
##' ## Joint effect of education and income
##' plot(fit1, which = c("education", "income"), n = 20)
##'
##' ## Bring up interactive menu
##' \dontrun{
##' plot(fit1, ask = TRUE)
##' 
##'   # displays menu:
##'   # Select one or two variable numbers (separated by spaces), or 0 to exit:
##' 
##'   # 1: education
##'   # 2: income
##'   # 3: type
##' }
plot.polywog <- function(x, which = NULL, ask = FALSE, auto.set.par = TRUE,
                         interval = TRUE, level = 0.95, bag = TRUE,
                         FUN3D = c("contour", "filled.contour", "wireframe",
                         "persp3d"),
                         control.plot = list(),
                         ...)
{
    FUN3D <- match.arg(FUN3D)
    pp <- NULL  # To avoid return error when no plot selected

    ## Extract regressor names
    xnames <- x$varNames

    if (ask) {
        while (TRUE) {
            ## Display menu of options
            askout <- paste(seq_along(xnames), ": ", xnames, sep="")
            askout <-
                c("Select one or two variable numbers (separated by spaces), or 0 to exit:",
                  "", askout, "")
            writeLines(askout)
            pick <- readline("Selection: ")

            ## Parse selection into variable names
            pick <- str_split(pick, " ")[[1]]
            pick <- as.integer(pick[pick != ""])
            if (any(is.na(pick)) || any(pick > length(xnames)))
                stop("Selected an invalid option")
            if (length(pick) > 2)
                stop("Cannot select more than two variables")
            if (any(pick < 1))
                break
            pick <- xnames[pick]
            pp <- preplotFromPick(x, pick = pick, interval = interval, level =
                                  level, bag = bag, ...)
            plot(pp, auto.set.par = auto.set.par, FUN3D = FUN3D, control.plot =
                 control.plot)
        }
    } else if (!is.null(which)) {
        if (length(which) > 2)
            stop("Cannot select more than two variables")
        pick <- if (is.numeric(which)) xnames[which] else which
        pp <- preplotFromPick(x, pick = pick, interval = interval, level =
                              level, bag = bag, ...)
        plot(pp, auto.set.par = auto.set.par, FUN3D = FUN3D, control.plot =
             control.plot)
    } else {
        ## Plot all univariate terms
        if (auto.set.par) {
            mfrow <- ceiling(sqrt(length(xnames)))
            mfcol <- ceiling(length(xnames) / mfrow)
            op <- par(mfrow = c(mfrow, mfcol))
            on.exit(par(op))
        }
        pp <- list()
        for (i in seq_along(xnames)) {
            ppi <- preplotFromPick(x, pick = xnames[i], interval = interval,
                                   level = level, bag = bag, ...)
            plot(ppi, auto.set.par = auto.set.par, FUN3D = FUN3D, control.plot =
                 control.plot)
            pp[[i]] <- ppi
        }
    }

    invisible(pp)
}