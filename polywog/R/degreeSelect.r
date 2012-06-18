##' @include helpers.r
##' @include polywog.r
NULL

degreeSelect <- function(...,
                         holdout.prop = 0.1,
                         max.degree = 5)
{
    cl <- match.call()

    ## Detect model formula
    if (is.null(cl[["formula"]])) {
        formula <- tryCatch(eval(cl[[2]], parent.frame()),
                            error = function(e) e)
        if (!inherits(formula, "formula") || inherits(formula, "tryError"))
            stop("'degreeSelect' must include a 'formula' argument")
    } else {
        formula <- eval(cl[["formula"]], parent.frame())
    }

    ## Assemble model frame
    formula <- as.Formula(formula)
    mf <- match(c("data", "subset", "weights", "na.action"), names(cl), 0L)
    mf <- cl[c(1L, mf)]
    mf$formula <- formula
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## Extract holdout sample
    nobs <- nrow(mf)
    ind.holdout <- sample(seq_len(nobs), ceiling(holdout.prop * nobs),
                          replace = FALSE)
    data.in <- mf[-ind.holdout, ]
    data.out <- mf[ind.holdout, ]

    degrees <- seq_len(max.degree)
    ans <- cbind(degree = degrees, "rms error" = NA)
    for (d in degrees) {
        ## Run polywog for given degree
        polywogCall <- cl
        polywogCall$holdout.prop <- polywogCall$max.degree <- NULL
        polywogCall$data <- data.in
        polywogCall$degree <- d
        polywogCall[[1]] <- as.name("polywog")
        polywogCall <- eval(polywogCall, parent.frame())

        ## Compute RMSE of out-of-sample fits
        pred <- predict(polywogCall, newdata = data.out)
        ans[d, 2] <- sqrt(mean((pred - model.response(data.out))^2))
    }

    ans
}

