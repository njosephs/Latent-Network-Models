#----------------------------------------
#
#               EM.LNM.U
#
#----------------------------------------

#pij -> dij 
logit <- function(p) log(p/(1 - p))
distance.u <- function(p) logit(1 - p/2)

#newton raphson 
nr.u <- function(x, y, s, nE, tol = 10e-4){
  #x is parameter to update 
  #y is other parameter in beta dist 
  # eta or phi
  #nE is number of edges 
  repeat{
    x.new <- x - ((digamma(x + y) -  digamma(x) +  sum(s)/nE) /(trigamma(x + y) - trigamma(x))) #nr update  
    if(abs((x.new - x)/x)<tol) break #convergence check
    x <- x.new #update new parameter est
  }
  return(x.new)
}

#EM
LNM.EM.U <- function(A, tol = 10e-4, no.iters = 1000){
  
  #An unweighted adj matrix
  
  #create Y vector -------------------------------------
  Y <- as.numeric(A[upper.tri(A)])
  nE <- length(Y) 
  
  #Initialize parameter values -------------------------
  alpha <- 1 #Prior beta values from uniform
  beta <- 1 #Prior beta values from uniform
  Q <- 0 #initialize Q
  iter <- 1
  
  #iterate until convergence
  repeat{
  
    #E Step 
    pi <- (alpha + Y)/(alpha + beta + 1) # p_ij
    eta <- digamma(alpha + Y) - digamma(alpha + beta+ 1) # log(p_ij)
    phi <- digamma(beta + 1 - Y) - digamma(alpha + beta + 1) # log(1-p_ij)
    
    #M Step 
    alpha <- nr.u(alpha, beta, eta, nE)
    beta <-  nr.u(beta, alpha, phi, nE)
    
    #Check convergence 
    Q.new <- sum((Y + alpha - 1)*eta  
                  - (Y - beta + 1)*phi 
                  - pi
                  +1 + log(gamma(alpha + beta)) - log(gamma(alpha)) - log(gamma(beta))
                  )
      
    if(iter > no.iters || (abs((Q.new - Q) / Q) < tol)) break
    Q <- Q.new
    iter <- iter + 1 
  }
  
  list(alpha = alpha, beta = beta, pi = pi, d = distance.u(pi), no.iter = iter)
}

#----------------------------------------
#
#               EM.LNM.W
#
#----------------------------------------

#pij -> dij 
distance.w <- function(l) 1/l 

#newton raphson 
nr.w <- function(alpha, beta, eta, nE, tol = 10e-4){
  #lapha is parameter to update 
  #beta, eta is other parameters in Gamma
  #nE is number of edges 
  repeat{
    alpha.new <- alpha - (digamma(alpha) - sum(eta)/nE - log(beta))/(trigamma(alpha)) #nr update  
    if(abs((alpha.new - alpha)/alpha)<tol) break #convergence check
    alpha <- alpha.new #update new parameter est
  }
  return(alpha.new)
}

#EM
LNM.EM.W <- function(W, tol = 10e-4, no.iters = 1000){
  
  #An unweighted adj matrix
  
  #create Y vector -------------------------------------
  Y <- as.numeric(W[upper.tri(W)])
  nE <- length(Y) 
  
  #Initialize parameter values -------------------------
  alpha <- 1 #Prior beta values from uniform
  beta <- 1 #Prior beta values from uniform
  Q <- 0 #initialize Q
  iter <- 1
  
  #iterate until convergence
  repeat{
    
    #E Step 
    pi <- (alpha + Y)/(beta + 1) # lambda_ij
    eta <- log(1 +beta) + digamma(alpha + Y) # log(lambda_ij)
    
    #M Step 
    beta <-  (nE * alpha) / (sum(pi))
    alpha <- nr.w(alpha, beta, eta, nE)
    
    
    #Check convergence 
    Q.new <- sum((Y + alpha - 1)*eta  
                 -(1 + beta)*pi 
                 - log(factorial(Y))
                 + alpha * log(beta) - log(gamma(alpha))
    )
    
    if(iter > no.iters || (abs((Q.new - Q) / Q) < tol)) break
    Q <- Q.new
    iter <- iter + 1 
  }
  list(alpha = alpha, beta = beta, pi = pi, d = distance.w(pi), no.iter = iter)
}





