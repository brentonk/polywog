## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model with bootstrap iterations
set.seed(22)
fit1 <- polywog(prestige ~ education + income + type, data = Prestige,
                boot = 10)

## Predicted prestige across occupational categories
predVals(fit1, "type")

## Predicted prestige by education
predVals(fit1, "education", n = 10)
