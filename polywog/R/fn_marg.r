################################################################################
### fn_marg.r
###
### Functions to compute pointwise and average marginal effects from polywog
### models.
###
### Brenton Kenkel and Curt Signorino
### created 2012-04-16
### last updated 2012-04-16 (bjk)
################################################################################

# NOT COMPLETE; DO NOT USE THESE FUNCTIONS YET

library(matrixStats)  # in package, should simply import 'rowProds'

## MAKE SURE TO TEST WITH MODELS WITH A SINGLE TERM INCLUDED LINEARLY, AND
## MULTIPLE TERMS INCLUDED LINEARLY

marginalEffects <- function(model, ...)
{
    ## Extract *raw* inputs (not polynomial-transformed) and matrix of which
    ## terms correspond to each coefficient
    X <- makeX(model$formula, model$model, 1)
    pt <- model$polyTerms

    ## Append linearly included terms, if any, to 'pt'
    ncoef <- length(coef(model)) - 1  # Remove intercept
    ndiff <- ncoef - nrow(pt)
    if (ndiff > 0) {
        pt <- rbind(pt, matrix(0, ndiff, ndiff))
        pt <- cbind(pt, rbind(matrix(0, nrow(pt), ndiff), diag(ndiff)))
    }

    ans <- matrix(0, nrow = nrow(X), ncol = ncol(X))

    ## Outer loop: each input variable
    for (i in seq_len(ncol(X))) {
        ## Inner loop: each coefficient
        for (j in seq_len(ncoef)) {
            ## Move on if not relevant to this variable
            if (pt[j, i] == 0)
                next

            ## Compute appropriate powers of input variables
            this_pt <- pt[j, ]
            this_pt[i] <- this_pt[i] - 1
            Xpow <- sweep(X, 2, this_pt, "^")

            ## Compute marginal effects (w.r.t. this coefficient)
            ans[, i] <- ans[, i] + (this_pt[i] + 1) * coef(model)[j] *
                rowProds(Xpow)
        }
    }
}
