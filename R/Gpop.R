Gpop.ttest.1sample <- function (xbar, sd, mu, n, tails = 1, alpha = 0.05, Gmax = 10, steps = 0.0001, full = F, digits = 5) {
  gate <- T
  c <- 1
  for (G in seq(from = 0, to = Gmax, by = steps)) {
    c = c + 1
  }
  results <- array(NA, c(c - 1, 3))

  c <- 1

  for (G in seq (from = 0, to = Gmax, by = steps)) {

      if (xbar > mu){
        mu_est <- mu + G * sd
        t  <- (xbar - mu_est) / (sd / sqrt(n))
        t0 <- (xbar - mu) / (sd / sqrt(n))
        p  <- tails*pt(-abs(t), df=n-1)
      } else {
        xbar_est <- xbar + G * sd
        t0 <- (xbar - mu) / (sd / sqrt(n))
        t  <- (xbar_est - mu ) / (sd / sqrt(n))
        p  <- tails*pt(-abs(t), df=n-1)
      }

    results [c, 1] = G
    results [c, 2] = t
    results [c, 3] = round(p, digits)

    c = c + 1
  }

  k <- 1
  for (k in seq (from = 1, to = c-2, by = 1)) {
    if (k < c - 2) {
      if (results [k, 3] >= 0.05) {
        if (results[k+1, 3] < 0.05) {
          if (gate == T){
            if (xbar > mu) {diff <- xbar - results[k, 1] * sd}
            else {diff <- xbar + results[k, 1] * sd}
            print (paste("P-value < ", alpha, " @ d = ", round(results[k, 1], digits), ", G =", round(results[k, 1] * sd, digits), ", xbar - G = ", round(diff, digits), ", p of G = ", round(tails*pt(-abs(results[k,2] - t0), df=n-1), digits)), quote = FALSE)
            if (full == F) {gate = F}
          }
        }
      }
      if (results [k, 3] <= 0.05) {
        if (results[k+1, 3] > 0.05) {
          if (gate == T){
            if (xbar > mu) {diff <- xbar - results[k, 1] * sd}
            else {diff <- xbar + results[k, 1] * sd}
            print (paste("P-value > ", alpha, " @ d = ", round(results[k, 1], digits), ", G =", round(results[k, 1] * sd, digits), ", xbar - G = ", round(diff, digits), ", p of G = ", round(tails*pt(-abs(results[k,2] - t0), df=n-1), digits)), quote = FALSE)
            if (full == F) {gate = F}
          }
        }
      }
    }
  }
  if (full == T) {
    return(results)
    assign("results", results, envir = .GlobalEnv)
  }
}

Gpop.ttest.2sample <- function (xbar, sd, n1, mu, sigma, n2, tails = 1, alpha = 0.05, Gmax = 10, steps = 0.0001, full = F, digits = 5) {
  gate <- T
  if (xbar < mu) {
    xbar0 <- mu
    mu0 <- xbar
    xbar <- xbar0
    mu <- mu0
    sd0 <- sigma
    sigma0 <- sd
    sd <- sd0
    sigma <- sigma0
  }
  c <- 1
  for (G in seq(from = 0, to = Gmax, by = steps)) {
    c = c + 1
  }
  results <- array(NA, c(c - 1, 3))

  c <- 1

  s1 = sd / sqrt(n1)
  s2 = sigma / sqrt(n2)

  for (G in seq (from = 0, to = Gmax, by = steps)) {

      if (xbar > mu){
        mu_est <- mu + G * sigma
        t  <- (xbar - mu_est) / sqrt(s1^2 / n1 + s2^2 / n2)
        t0 <- (xbar - mu) / sqrt(s1^2 / n1 + s2^2 / n2)
        p  <- tails*pt(-abs(t), df=n1 + n2 - 2)

      }
      else{
        xbar_est <- xbar - G * sd
        t  <- (xbar_est - mu ) / sqrt(s1^2 / n1 + s2^2 / n2)
        t0 <- (xbar - mu) / sqrt(s1^2 / n1 + s2^2 / n2)
        p  <- tails*pt(-abs(t), df=n1 + n2 - 2)
      }


    results [c, 1] = G
    results [c, 2] = t
    results [c, 3] = round(p,digits)

    c = c + 1

  }

  k <- 1
  for (k in seq (from = 1, to = c-2, by = 1)) {
    if (k < c - 2) {
      if (results [k, 3] >= 0.05) {
        if (results[k+1, 3] < 0.05) {
          if (gate == T){
            if (xbar > mu) {
              diff <- round(xbar - mu + results[k, 1] * sigma, digits)
              Res_G <- round(results[k, 1] * sigma, digits)
              }
            print (paste("P-value > ", alpha, " @ d = ", round(results[k, 1],digits), ", G =", Res_G, ", xbar - mu + G = ", diff, ", p of G = ", round(tails*pt(-abs(results[k,2] - t0), df=n1-1),digits)), quote = FALSE)
            if (full == F) {gate = F}
            }
        }
      }
      if (results [k, 3] <= 0.05) {
        if (results[k+1, 3] > 0.05) {
          if (gate == T){
            if (xbar > mu) {
              diff <- round(xbar - mu + results[k, 1] * sigma, digits)
              Res_G <- round(results[k, 1] * sigma, digits)
            }
            if (xbar < mu) {

              diff <- round(xbar - mu + results[k, 1] * sd, digits)
              Res_G <- round(results[k, 1] * sd, digits)
            }
            print (paste("P-value > ", alpha, " @ d = ", round(results[k, 1],digits), ", G =", Res_G, ", xbar - mu + G = ", diff, ", p of G = ", round(tails*pt(-abs(results[k,2] - t0), df=n1-1),digits)), quote = FALSE)
            if (full == F) {gate = F}
          }
        }
      }
    }
  }
  if (full == T) {
    return(results)
    assign("results", results, envir = .GlobalEnv)
  }
}

Gpop.regression <- function (B, SE, n, k, tails = 1, alpha = 0.05, Gmax = 10, steps = 0.0001, full = F, digits = 5) {
  gate <- T
  c <- 1
  for (G in seq(from = 0, to = Gmax, by = steps)) {
    c = c + 1
  }
  results <- array(NA, c(c - 1, 3))

  c <- 1

  for (G in seq (from = 0, to = Gmax, by = steps)) {

    mu_est <- abs(G * B)
    t  <- (abs(B) - mu_est) / SE
    t0 <- abs(B) / SE
    p  <- tails*pt(-abs(t), df=n-k-1)

    results [c, 1] = G
    results [c, 2] = t
    results [c, 3] = round(p,digits)

    c = c + 1

  }

  q <- 1

  for (q in seq (from = 1, to = c-2, by = 1)) {
    if (q < c - 2) {
      if (results [q, 3] >= 0.05) {
        if (results[q+1, 3] < 0.05) {
          if (gate == T){
            diff <- round(abs(B) - abs(results[q, 1] * B), digits)
            Res_c <-round(results[q, 1], digits)
            Res_G <- round(abs(results[q, 1] * B), digits)
            Frank <- abs(round((1 - results[q, 1] * B)/B, digits))
            print (paste("P-value < ", alpha, " @ c = ", Res_c, ", G =", Res_G, ", B - G = ", diff, ", p of G = ", round(tails*pt(-abs(results[q,2] - t0), df=n-k-1), digits)), quote = FALSE)
          }
        }
      }
      if (results [q, 3] <= 0.05) {
        if (results[q+1, 3] > 0.05) {
          diff <- round(abs(B) - abs(results[q, 1] * B), digits)
          Res_c <-round(results[q, 1], digits)
          Res_G <- round(abs(results[q, 1] * B), digits)
          Frank <- abs(round((1 - results[q, 1] * B)/B, digits))
          print (paste("P-value > ", alpha, " @ c = ", Res_c, ", G =", Res_G, ", B - G = ", diff, ", p of G = ", round(tails*pt(-abs(results[q,2] - t0), df=n-k-1), digits)), quote = FALSE)
          if (full == F) {gate = F}
        }
      }
    }
  }
  if (full == T) {
    return (results)
    assign("results", results, envir = .GlobalEnv)
  }
}
