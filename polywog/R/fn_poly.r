################################################################################
### fn_poly.r
###
### Functions to create polynomial expansions of matrices
###
### Brenton Kenkel and Curt Signorino
### created 2011-10-25
### last updated 2012-04-30 (bjk)
################################################################################

repAsMat <- function(x, n = 1)
{
    matrix(x, nrow = n, ncol = length(x), byrow = TRUE)
}

polyCombos <- function(k, degree)
{
    y <- diag(k)
    if (degree < 2)
        return(data.frame(y))
    ans <- vector("list", degree)
    ans[[1]] <- y

    for (i in 2:degree) {
        y <- vector("list", nrow(ans[[i-1]]))
        for (j in 1:nrow(ans[[i-1]]))
            y[[j]] <- repAsMat(ans[[i-1]][j, ], k) + diag(k)
        y <- do.call("rbind", y)
        ans[[i]] <- y
    }

    ans <- data.frame(do.call("rbind", ans))
    ans <- ans[!duplicated(ans), ]
    return(ans)
}

##
## based on function "poly" in R's stats package, used under GPL
##
## modified version used to avoid warnings about degree being less than number
## of unique points
## 

##' Calculate raw polynomials
##'
##' This is essentially the same as running \code{poly(x, degree = degree, raw =
##' TRUE)} (and is based on the \code{\link{poly}} code), but \code{NA}s are
##' allowed and no warning is issued when the degree exceeds the number of
##' unique points.
##' @param x a numeric vector.
##' @param degree the degree of the polynomial.
##' @param na.ok logical: whether to allow \code{NA}s.  If \code{na.ok = FALSE},
##' then the function will stop with an error if any \code{NA}s are present in
##' \code{x}.
##' @return A matrix containing a raw polynomial expansion of \code{x}.
##' @author Brenton Kenkel and Curtis S. Signorino
##' @export
##' @examples
##' x <- c(1, 0, 1, 0, 0, NA, 1)
##' rawpoly(x, degree = 3, na.ok = TRUE)
rawpoly <- function(x, degree = 1, na.ok = FALSE)
{
    if (degree < 1)
        stop("'degree' must be at least 1")
    if (!na.ok && any(is.na(x)))
        stop("missing values are not allowed if 'na.ok' is FALSE")

    Z <- outer(x, 1L:degree, "^")
    colnames(Z) <- 1L:degree
    attr(Z, "degree") <- 1L:degree
    class(Z) <- c("poly", "matrix")
    
    return(Z)
}

##' Calculate polynomial basis expansion of a matrix
##'
##' This produces the same output as running \code{polym(x, degree = degree, raw
##' = TRUE)} (and is based on the \code{\link{polym}} code).  The differences
##' are: \itemize{
##'   \item Computation is much faster than \code{polym} when \code{degree} is
##' relatively low and \code{x} contains many columns (the typical use case for
##' \code{\link{polywog}}).
##'   \item Column names of the output use the column names of \code{x} in order
##' to be more intelligible.
##'   \item \code{NA}s are allowed via the \code{na.ok} option.
##' }
##' @param x a numeric matrix.
##' @inheritParams rawpoly
##' @return A numeric matrix containing the specified basis expansion of the
##' input matrix, along with attributes \code{degree} giving the degree of each
##' column and \code{polyTerms} giving the power of each input term represented
##' in each column.
##' @author Brenton Kenkel
##' @export
##' @examples
##' x <- matrix(rnorm(20), ncol = 2)
##' colnames(x) <- c("a", "b")
##' polym2(x, degree = 3)
polym2 <- function (x, degree = 1, na.ok = FALSE)
{
    x <- unclass(as.data.frame(x))

    ## This block is taken straight from the stats::polym source
    nd <- length(x)
    if (nd == 0) 
        stop("must supply one or more vectors")
    if (nd == 1) 
        return(rawpoly(x[[1L]], degree, na.ok))
    n <- sapply(x, length)
    if (any(n != n[1L])) 
        stop("arguments must have the same length")

    ## New method for generating possible combinations, as the default was
    ## ridiculously inefficient for matrices with many columns
    z <- polyCombos(nd, degree)
    s <- rowSums(z)
    
    ## This block is also straight from stats::polym
    res <- cbind(1, rawpoly(x[[1L]], degree, na.ok))[, 1 + z[, 1]]
    for (i in 2:nd)
        res <- res * cbind(1, rawpoly(x[[i]], degree, na.ok))[, 1 + z[, i]]

    ## Making intelligible column names
    names(z) <- names(x)
    oldz <- z
    for (i in seq_len(nd)) {
        zz <- z[, i]
        z[, i] <- as.character(z[, i])
        z[, i] <- paste(names(x)[i], "^", z[, i], sep = "")
        z[zz == 0, i] <- NA
        z[zz == 1, i] <- names(x)[i]
    }
    colnames(res) <- apply(z, 1L,
                           function(x) paste(x[!is.na(x)], collapse = "."))
    attr(res, "degree") <- as.vector(s)
    attr(res, "polyTerms") <- as.matrix(oldz, rownames.force = FALSE)
    return(res)
}
