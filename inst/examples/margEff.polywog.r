## Using occupational prestige data
data(Prestige, package = "car")
Prestige <- transform(Prestige, income = income / 1000)

## Fit a polywog model
set.seed(22)
fit1 <- polywog(prestige ~ education + income | type, data = Prestige)

## Compute marginal effects for all variables
me1 <- margEff(fit1)
summary(me1)  # type was included linearly, hence constant effects

## Plotting density of the results
plot(me1)

## Can do the same when just examining a single variable
me2 <- margEff(fit1, xvar = "income")
summary(me2)
plot(me2)
