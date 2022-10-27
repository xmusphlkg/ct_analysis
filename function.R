fit_best <- function(data, distribution.type = NULL){
  data[data == 0] <- 1/2
  fit.gamma <- try(fitdistr(data, "gamma"))
  if (class(fit.gamma) == "try-error") {
    fit.gamma$loglik <- NA
  }
  fit.weib <- try(fitdistr(data, "weibull"))
  if (class(fit.weib) == "try-error") {
    fit.weib$loglik <- NA
  }
  fit.lognorm <- try(fitdistr(data, "log-normal"))
  if (class(fit.lognorm) == "try-error") {
    fit.lognorm$loglik <- NA
  }
  fit.type <- c("gamma", "weibull", "lognormal")
  distribution.type <- ifelse(is.null(distribution.type),
                              fit.type[which.max(c(fit.gamma$loglik,
                                                   fit.weib$loglik, fit.lognorm$loglik))],
                              distribution.type)
  x <- NULL
  rm(x)
  if (distribution.type == "gamma") {
    shape <- fit.gamma$estimate[1]
    rate <- fit.gamma$estimate[2]
    mean <- shape/rate
    sd <- sqrt(shape)/rate
    return(list(distr = 'gamma',
                shape = shape,
                rate = rate,
                mean = mean,
                sd = sd))
  }
  else if (distribution.type == "weibull") {
    shape <- fit.weib$estimate[1]
    scale <- fit.weib$estimate[2]
    mean <- scale * exp(lgamma(1 + 1/shape))
    sd <- sqrt(scale^2 * (exp(lgamma(1 + 2/shape)) - (exp(lgamma(1 + 1/shape)))^2))
    return(list(distr = 'weibull',
                shape = shape,
                scale = scale,
                mean = mean,
                sd = sd))
  }
  else if (distribution.type == "lognormal") {
    meanlog <- fit.lognorm$estimate[1]
    sdlog <- fit.lognorm$estimate[2]
    mean <- exp(1/2 * sdlog^2 + meanlog)
    sd <- sqrt(exp(2 * meanlog + sdlog^2) * (exp(sdlog^2) - 1))
    return(list(distr = 'lognormal',
                meanlog = meanlog,
                sdlog = sdlog,
                mean = mean,
                sd = sd))
  }
}

expose_date <- function(a, b){
  return(as.Date(seq.Date(a, b, by = 'day')))
}

expose_seq <- function(x, data){
  infector_expose_dates <- seq.Date(from = data[x, 'infector_exposedate1'],
                                    to = data[x, 'infector_exposedate2'],
                                    by = 'day')
  infectee_expose_dates <- seq.Date(from = data[x, 'infectee_exposedate1'],
                                    to = data[x, 'infectee_exposedate2'],
                                    by = 'day')
  datafile <- as.data.frame(expand.grid(infectee_expose_dates, infector_expose_dates)) %>%
    mutate(median = as.numeric(Var1 - Var2),
           median = ifelse(median<0, NA, median))
  gt_median <- median(datafile$median, na.rm = T)
  gt_max <- max(infectee_expose_dates) - min(infector_expose_dates)
  gt_min <- min(infectee_expose_dates) - max(infector_expose_dates)
  gt_min <- ifelse(gt_min<0, 0, gt_min)
  return(c(gt_min, gt_median, gt_max))
}

fit_best <- function(data, distribution.type = NULL){
  data[data == 0] <- 1/2
  fit.gamma <- try(fitdistr(data, "gamma"))
  if (class(fit.gamma) == "try-error") {
    fit.gamma$loglik <- NA
  }
  fit.weib <- try(fitdistr(data, "weibull"))
  if (class(fit.weib) == "try-error") {
    fit.weib$loglik <- NA
  }
  fit.lognorm <- try(fitdistr(data, "log-normal"))
  if (class(fit.lognorm) == "try-error") {
    fit.lognorm$loglik <- NA
  }
  fit.type <- c("gamma", "weibull", "lognormal")
  distribution.type <- ifelse(is.null(distribution.type),
                              fit.type[which.max(c(fit.gamma$loglik,
                                                   fit.weib$loglik, fit.lognorm$loglik))],
                              distribution.type)
  x <- NULL
  rm(x)
  if (distribution.type == "gamma") {
    shape <- fit.gamma$estimate[1]
    rate <- fit.gamma$estimate[2]
    mean <- shape/rate
    sd <- sqrt(shape)/rate
    return(list(distr = 'gamma',
                shape = shape,
                rate = rate,
                mean = mean,
                sd = sd))
  }
  else if (distribution.type == "weibull") {
    shape <- fit.weib$estimate[1]
    scale <- fit.weib$estimate[2]
    mean <- scale * exp(lgamma(1 + 1/shape))
    sd <- sqrt(scale^2 * (exp(lgamma(1 + 2/shape)) - (exp(lgamma(1 + 1/shape)))^2))
    return(list(distr = 'weibull',
                shape = shape,
                scale = scale,
                mean = mean,
                sd = sd))
  }
  else if (distribution.type == "lognormal") {
    meanlog <- fit.lognorm$estimate[1]
    sdlog <- fit.lognorm$estimate[2]
    mean <- exp(1/2 * sdlog^2 + meanlog)
    sd <- sqrt(exp(2 * meanlog + sdlog^2) * (exp(sdlog^2) - 1))
    return(list(distr = 'lognormal',
                meanlog = meanlog,
                sdlog = sdlog,
                mean = mean,
                sd = sd))
  }
}

