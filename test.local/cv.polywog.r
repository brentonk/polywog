## Test of the cross-validation function

## Generate a model of degree 2
n <- 1000
k <- 8
X <- matrix(rnorm(n * k), ncol = k)
colnames(X) <- head(letters, k)
pt <- makePolyTerms(2, k - 4, 4, rep(FALSE, k), colnames(X))
X.expand <- expandMatrix(X, pt)
b <- rnorm(nrow(pt)) * rbinom(nrow(pt), 1, 1/2)
names(b) <- rownames(pt)
y <- 1 + as.numeric(X.expand %*% b) + rnorm(n)
dat <- data.frame(y, X)

## Run cv.polywog with degrees 1-4
cv_alasso <- cv.polywog(y ~ ., data = dat, degrees.cv = 1:4)

## Do the same when fitting via SCAD
cv_scad <- cv.polywog(y ~ ., data = dat, method = "scad", degrees.cv = 1:4)
