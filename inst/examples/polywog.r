## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model with bootstrap iterations
set.seed(22)
fit1 <- polywog(prestige ~ education + income + type, data = Prestige,
                boot = 10)

## Basic information
print(fit1)
summary(fit1)

## See how fitted values change with education holding all else fixed
predVals(fit1, "education", n = 10)

## Plot univariate relationships
plot(fit1)

## Change regularization method
fit2 <- update(fit1, method = "scad")
cbind(coef(fit1), coef(fit2))
