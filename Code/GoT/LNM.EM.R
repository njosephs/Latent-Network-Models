#----------------------------------------
#
#               EM.LNM
#
#----------------------------------------

#pij -> dij 
logit <- function(p) log(p/(1 - p))
distance <- function(p) logit(1 - p/2)

#newton raphson 
nr <- function(x, y,theta, nE, tol = 10e-4){
  #x is parameter to update 
  #y is other parameter in beta dist 
  #theta weighted sum of either eta or gamma
  #nE is number of edges 
  repeat{
    x.new <- x - ((digamma(x + y) -  digamma(x) +  sum(eta)/nE) /(trigamma(x + y) - trigamma(x))) #nr update  
    if(abs((x.new - x)/x)<tol) break #convergence check
    x <- x.new #update new parameter est
  }
  return(x.new)
}

#EM
LNM.EM <- function(A, tol = 10e-4, no.iters = 1000){
  
  #create Y vector -------------------------------------
  Y <- as.numeric(A[upper.tri(A)])
  nE <- length(Y) 
  
  #Initialize parameter values -------------------------
  alpha <- 1 #Prior beta value
  beta <- 1 #Prior beta value
  Q <- 0 #initialize 
  
  #iterate until convergence
  repeat{
  
    #E Step 
    pi <- (alpha + Y)/(alpha + beta + 1) # p_ij
    eta <- digamma(alpha + Y) - digamma(alpha + beta+ 1) # log(p_ij)
    gamma <- digamma(beta + 1 - Y) - digamma(alpha + beta + 1) # log(1-p_ij)
    
    #M Step 
    alpha.new <- nr(alpha, beta, eta, nE)
    beta.new <-  nr(beta, alpha, gamma, nE)
    
    #Check convergence 
    Q.new <- sum((Y + alpha.new - 1)*eta  
                  - (Y - beta.new + 1)*gamma 
                  - pi
                  +1 + log(gamma(alpha.new + beta.new)) - log(gamma(alpha.new)) - log(gamma(beta.new))
                  )
      
    if (abs((Q.new - Q) / Q) < tol) break
    Q <- Q.new
  }
  
  list(alpha = alpha, beta = beta, pi = pi, d = distance(pi))
}




