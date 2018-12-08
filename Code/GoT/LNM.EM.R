#----------------------------------------
#
#               EM.LNM
#
#----------------------------------------

#pij -> dij 
logit <- function(p) log(p/(1 - p))
distance <- function(p) logit(1 - p/2)

#newton raphson 
nr <- function(x, y, s, nE, tol = 10e-4){
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
LNM.EM <- function(A, tol = 10e-4, no.iters = 1000){
  
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
    alpha <- nr(alpha, beta, eta, nE)
    beta <-  nr(beta, alpha, phi, nE)
    
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
  
  list(alpha = alpha, beta = beta, pi = pi, d = distance(pi), no.iter = iter)
}



