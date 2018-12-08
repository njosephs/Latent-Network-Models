LOGEPS <- log(.Machine$double.eps / 2)

lse <- function (x) {
  m <- max(x); x <- x - m
  m + log(sum(exp(x[x > LOGEPS])))
}

lse2 <- function (x, y) {
  m <- pmax(x, y); d <- -abs(x - y)
  ifelse(d < LOGEPS, m, m + log(1 + exp(d))) }

rlcat <- function (n, l) {
  l <- Reduce(lse2, l, accumulate = TRUE) # "cumlse" 
  l <- l - l[length(l)] # normalize 
  findInterval(log(runif(n)), l) + 1
}

LNM.MCMC <- function(G, Nk = 2, d = 2, ns = 10000) {
  library(MASS)
  library(invgamma)
  
  Nv <- nrow(G)
  
  # Initialize
  mu <- array(dim = c(Nk, d, ns)); mu[, , 1] <- 0
  sigma <- matrix(nrow = Nk, ncol = ns); sigma[, 1] <- 1
  lambda <- matrix(nrow = Nk, ncol = ns); lambda[, 1] <- -log(Nk)
  K <- matrix(nrow = Nv, ncol = ns); K[, 1] <- sample.int(Nk, Nv, replace = TRUE)
  Z <- array(dim = c(Nv, d, ns)); Z[, , 1] <- 0
  
  # Updates
  for (t in 2:ns) {
    
    # Mu
    for (k in 1:Nk) {
      Zg <- apply(matrix(Z[K[, t - 1] == k, , t - 1], ncol = d), 2, mean)
      Zg <- ifelse(is.na(Zg), 0, Zg) # in case no one is in group k
      ng <- sum(K[, t - 1] == k)
      mu[k, , t] <- mvrnorm(n = 1
                            , mu = ( ng * Zg ) / ( ng + sigma[k, t - 1] ) 
                            , Sigma = ( (sigma[k, t -1 ]) /  (ng + sigma[k, t -1]) ) * diag(d))
    }
    
    # Sigma
    for (k in 1:Nk) {
      ng <- sum(K[, t - 1] == k) 
      SS_Zg <- sum(apply(matrix(Z[K[, t - 1] == k, , t - 1] - mu[k, , t], ncol = d), 2, crossprod))
      sigma[k, t] <- (1 + SS_Zg) * rinvchisq(1, 1 + ng*d)
    }
    
    # Group K
    for (i in 1:Nv) {
      lambda_g <- sapply(1:Nk, FUN = function(k) {
        C <- chol(sigma[k, t] * diag(d))
        y <- backsolve(C, (Z[i, , t - 1] - mu[k, , t]), transpose = TRUE)
        - (1/2) * log(2 * pi) - sum(log(diag(C))) - sum(y^2) / 2
      })
      
      lambda[, t] <- lambda_g - lse(lambda_g)
      K[i, t] <- rlcat(1, lambda[, t]) 
    }
    
    # Latent variable Z
    Z[, , t] <- Z[, , t - 1]
    for (i in 1:Nv) {
      
      Zstar <- mvrnorm(1, mu = rep(0, d), Sigma = diag(d))
      
      C <- chol(sigma[K[i, t - 1], t] * diag(d))
      y <- backsolve(C, (Z[i, , t - 1] - mu[K[i, t], , t]), transpose = TRUE)
      y.star <- backsolve(C, (Zstar - mu[K[i, t], , t]), transpose = TRUE)
      logR <- (-(1/2) * log(2 * pi) - sum(log(diag(C)))- sum(y^2) / 2) - # p(Z*)
        (-(1/2) * log(2 * pi) - sum(log(diag(C))) - sum(y^2) / 2)        # p(Z^(t-1))
      
      if (logR >= 0 || log(runif(1)) < logR) Z[i, , t] <- Zstar
    }
  }
  return(list(mu = mu, sigma = sigma, lambda = lambda, K = K, Z = Z)) 
}
  