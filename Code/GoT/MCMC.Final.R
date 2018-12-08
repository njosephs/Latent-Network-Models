LNM.MCMC <- function(G, Nk = 2, d = 2, ns = 10000) {
  Nv <- nrow(G)
  
  # Initialize
  mu <- array(dim = c(Nk, d, ns)); mu[, , 1] <- 0
  sigma <- matrix(nrow = Nk, ncol = ns); sigma[, 1] <- 1
  lambda <- matrix(nrow = Nk, ncol = ns); lambda[, 1] <- 1/Nk
  K <- matrix(nrow = Nv, ncol = ns); K[, 1] <- sample.int(Nk, Nv, replace = TRUE)
  Z <- array(dim = c(Nv, d, ns)); Z[, , 1] <- 0
  
  # Updates
  for (t in 2:ns) {
    
    # Group K
    K[, t] <- sample.int(Nk, Nv, replace = TRUE)
    
    # Mu
    for (k in 1:Nk) {
      mu[k, , t] <- mvrnorm(n = 1
                            , mu = Z[k, , t]/(sigma[k, ] + 1)
                            , Sigma = (sigma[k, ] / (sigma[k, ] + diag(x = 1, nrow = 2, ncol = 2))) )
    }
    
    # Sigma
    for (k in 1:Nk) {
      sigma[k, t] # DO
    }
    
    # Latent variable Z
    Z[, t] <- Z[, t - 1]
    for (i in 1:Nv) {
      
      Zstar # DO
      logR # DO
      
      if (logR > 0 || log(runif(1)) < logR) Z[i, t] <- Zstar
    }
    
    return() # DO
}