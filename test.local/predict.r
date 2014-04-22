## Simulate data and coefficient vector
n <- 1000
k <- 10
X <- matrix(rnorm(n * k), ncol = k)
pt <- makePolyTerms(degree = 2, k_expand = k, k_lin = 0,
                    binary_cols = rep(FALSE, k), names. = letters[1:k])
Xe <- expandMatrix(X, pt, TRUE)
b <- rnorm(ncol(Xe))

## Predicted values
pv_man <- as.numeric(Xe %*% b)
pv_cpp <- computePredict(X = X,
                         poly_terms = pt,
                         coef = list(main = b),
                         forPredVals = FALSE,
                         interval = FALSE,
                         bag = FALSE,
                         level = 0,
                         logit = FALSE)$fit
all.equal(pv_man, pv_cpp)

## Predicted probabilities
pp_man <- plogis(pv_man)
pp_cpp <- computePredict(X = X,
                         poly_terms = pt,
                         coef = list(main = b),
                         forPredVals = FALSE,
                         interval = FALSE,
                         bag = FALSE,
                         level = 0,
                         logit = TRUE)$fit
all.equal(pp_man, pp_cpp)
