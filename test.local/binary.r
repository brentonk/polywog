## Test that everything works with binary models

n <- 10000
X <- matrix(rnorm(4 * n), ncol = 4)
colnames(X) <- letters[1:4]
ystar <- 1 + X[, 1] - X[, 2]^2 + X[, 1] * X[, 3]^2 + rlogis(n)
y <- as.numeric(ystar > 0)
dat <- data.frame(y, X)

## Logistic model via adaptive LASSO
fit_alasso_lm <- polywog(y ~ ., data = dat, family = "binomial")
fit_alasso_lm <- bootPolywog(fit_alasso_lm, nboot = 5, thresh = 1e-4)
fit_alasso_glm <- update(fit_alasso_lm, penwt.method = "glm")
fit_alasso_glm <- bootPolywog(fit_alasso_glm, nboot = 5, thresh = 1e-4)

## Logistic model via SCAD
fit_scad <- update(fit_alasso_lm, method = "scad")
fit_scad <- bootPolywog(fit_scad, nboot = 5)

## Predicted probabilities
pred_link <- predict(fit_alasso_lm, type = "link", interval = TRUE)
pred_response <- predict(fit_alasso_lm, type = "response", interval = TRUE)
all(plogis(pred_link[, "fit"]) == pred_response[, "fit"])
