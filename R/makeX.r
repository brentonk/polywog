##
## Assemble the 1-dimensional (i.e., non-polynomial-expanded) model matrix
##
## ARGUMENTS:
##   formula: model formula (of class "Formula")
##   mf: model frame
##
## RETURN:
##   model matrix (no intercept)
##
makePlainX <- function(formula, mf)
{
    ## Ensure no intercepts included
    formula <- removeIntercepts(formula, mf)

    ## Extract the model matrix from the model frame
    X <- model.matrix(formula, data = mf, rhs = 1)
    if (length(formula)[2] > 1)
        X <- cbind(X, model.matrix(formula, data = mf, rhs = 2))

    X
}

##
## Assemble polynomial-expanded model matrix (plus linear terms if specified)
## without intercept.
##
## This function is used within 'polywog' and some post-estimation analysis
## functions; it typically should not be called directly by a user.
##
## ARGUMENTS:
##   formula: model formula (of class "Formula")
##   mf: model frame
##   degree: degree of polynomial expansion
##   na.ok: whether to let NAs pass through
##
## RETURN:
##   model matrix
##
makeX <- function(formula, mf, degree, na.ok = FALSE)
{
    ## Ensure no intercepts included
    formula <- removeIntercepts(formula, mf)

    ## Polynomial expansion
    X <- model.matrix(formula, data = mf, rhs = 1)
    X <- polym2(X, degree = degree, na.ok = na.ok)

    ## Linear terms (if any)
    if (length(formula)[2] > 1) {
        ## Need to save "polyTerms" attribute or cbind() will destroy it
        pt1 <- attr(X, "polyTerms")
        Xlin <- model.matrix(formula, data = mf, rhs = 2)
        X <- cbind(X, Xlin)

        ## Add the new variables to "polyTerms"
        nlin <- ncol(Xlin)
        pt <- cbind(pt1, matrix(0, nrow = nrow(pt1), ncol = nlin))
        pt <- rbind(pt,
                    cbind(matrix(0, nrow = nlin, ncol = ncol(pt1)), diag(nlin)))
        colnames(pt) <- c(colnames(pt1), colnames(Xlin))
        attr(X, "polyTerms") <- pt
    }

    return(X)
}
