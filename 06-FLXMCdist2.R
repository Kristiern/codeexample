# Driver to be used with R package flexmix
# Written using the driver FLXMCdist1 as base
# For master's thesis
# Kristi Ernits

FLXMCdist2 <- function(formula=.~., dist, ...) {
  foo <- paste("FLXMC", dist, sep ="")
  if (!exists(foo))
    stop("This distribution has not been implemented yet.")
  get(foo)(formula, ...)
}

prepoc.y.pos.1 <- function(x) {
  if (ncol(x) > 1)
    stop("y must be univariate")
  if (any(x <= 0))
    stop("values must be > 0")
  x
}

FLXMCrayleigh <- function(formula=.~., ...)
{
  z <- new("FLXMC", weighted = TRUE, formula = formula,
           dist = "rayleigh", name = "model-based Rayleigh clustering")
  z@preproc.y <- prepoc.y.pos.1
  
  z@defineComponent <- function(para) {
    predict <-  function(x, ...)
      matrix(para$scale, nrow = nrow(x), ncol = 1, byrow = TRUE)
    
    logLik <- function(x, y)
      VGAM::drayleigh(y, scale = predict(x, ...), log = TRUE)
    
    new("FLXcomponent", parameters = list(scale = para$scale),
        predict = predict, logLik = logLik, df = para$df)
  }
  
  z@fit <- function(x, y, w, ...) {
    meanw <- mean(w)
    scale <- sqrt(mean(w * (y)**2) / (2*meanw))
    z@defineComponent(list(scale = scale, df = 1))
  }
  z
}



FLXMCtnorm0 <- function(formula=.~., method = "Nelder-Mead", warn = -1, ...)
{
  z <- new("FLXMC", weighted = TRUE, formula = formula,
           name = "model-based truncated normal clustering",
           dist = "tnorm0")
  z@preproc.y <- prepoc.y.pos.1
  
  z@defineComponent <- function(para) {
    predict <- function(x, ...) 
      matrix(para$mu, nrow = nrow(x), ncol = length(para$mu),
             byrow = TRUE)
    
    logLik <- function(x, y, ...) 
      log(truncnorm::dtruncnorm(y, a=0, b=Inf, mean = predict(x, ...), sd = para$scale))
    
    new("FLXcomponent", parameters = list(mu = para$mu, scale = para$scale),
        predict = predict, logLik = logLik, df = para$df)
  }
  
  z@fit <- function(x, y, w, component){
    if (!length(component)) {
      sw <- sum(w)
      mean <- sum(y * w) / sw
      var <- (sum(y^2 * w) / sw - mean^2) * sw / (sw - 1)
      mu <- mean
      scale <- sqrt(var)
      start <- c(mu, scale)
    } else start <- unname(unlist(component))
    f <- function(parms) -sum(log(truncnorm::dtruncnorm(y, a=0, b=Inf, mean = parms[1], sd = parms[2])) * w)
    oop <- options(warn = warn)
    on.exit(oop)
    parms <- optim(start, f, method = method)$par
    z@defineComponent(list(mu = parms[1], scale = parms[2], df = 2))
  }
  z
}

