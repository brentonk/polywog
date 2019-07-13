
logLik.polywog <- function(object, ...){
  if(is.null(object$y)){
    object <- update(object, y=TRUE)
  }
  y <- object$y
  eta <- object$fitted
  yhat <- do.call(object$family, list())$linkinv(eta)
  if(object$family == "gaussian"){
    p <- sum(object$coefficients != 0)
    e <- y - yhat
    ll <- dnorm(yhat, y, sqrt(sum(e^2)/(length(y))), log=TRUE)
    out <- sum(ll)
      }
  if(object$family == "binomial"){
    ll <- dbinom(y, 1, yhat, log=TRUE)
    out <- sum(ll)
  }
  if(!(object$family %in% c("gaussian", "binomial"))){
    stop("Currently, logLik.polywog only works for gaussian and binomial families.\n")
  }  
  return(out)
}

AIC.polywog <- function(object){
  p <- sum(object$coefficients != 0)
  ll <- logLik(object)
  -2*ll + 2*p
}
BIC.polywog <- function(object){
  p <- sum(object$coefficients != 0)
  ll <- logLik(object)
  -2*ll + p*log(length(object$fitted))
}
