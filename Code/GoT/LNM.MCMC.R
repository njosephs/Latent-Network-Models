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
      z.bar.g <- apply(Z[K[, t-1] == k, , t - 1], 2, mean)
      ng <- sum(K[, t-1] == k)
      mu[k, , t] <- mvrnorm(n = 1
                            , mu = ( ng * z.bar.g ) / ( ng + sigma[k, t -1] ) 
                            , Sigma = ((sigma[k, t -1 ]) /  (ng + sigma[k, t -1])) * diag(d))
    }
    
    # Sigma
    for (k in 1:Nk) {
      ng <- sum(K[, t-1] == k) 
      SS_zg <- sum(apply(Z[K[, t-1] == k, , t -1] - mu[k, , t], 2, crossprod))
      sigma[k, t] <- (1 + SS_zg) * rinvchisq(1, 1 + ng*d)
    }
    
    # Group K
    l_p <- sapply(1: Nk, FUN = function(k) {
    C <- chol(Sigma[k,t] * diag(d))
    y <- backsolve(C, (Z[i, , t-1] - mu[k, ,t]), transpose = TRUE)
    log_P <- - (1/2) * log(2* pi) - sum(log(diag(C)))- sum(y ^2) /2
    })
    
    lambda[, t] <- l_p - lse(l_p)
    K[, t] <- rlcat(Nv, lambda[, t])
    
    # Latent variable Z
    Z[, t] <- Z[, t - 1]
    for (i in 1:Nv) {
      
      Zstar <- mvrnorm(rep(0, d), diag(d))
      
      C <- chol(Sigma[K[i, t-1],t] * diag(d))
      y <- backsolve(C, (Z[i, , t-1] - mu[K[i, t], ,t]), transpose = TRUE)
      log_P <- - (1/2) * log(2* pi) - sum(log(diag(C)))- sum(y ^2) /2
    
      y.star <- backsolve(C, (Zstar[i, , t-1] - mu[K[i, t], ,t]), transpose = TRUE)
      log_P.star <- - (1/2) * log(2* pi) - sum(log(diag(C)))- sum(y ^2) /2
      logR <- log_P.star - log_P
      
      if (logR > 0 || log(runif(1)) < logR) Z[i, t] <- Zstar
    }
  }
  return(list(mu = mu, sigma = sigma, lambda = lambda, K = K, Z = Z)) 
}
  