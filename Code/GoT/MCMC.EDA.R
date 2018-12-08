load("./Desktop/GitHub/Latent-Network-Models/Data/A.Rdata")
source("./Desktop/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

ns <- 10000
Nv <- nrow(A)
Nk <- 2
d <- 2
nchains <- 1
params <- c(paste0("mu", 1:Nk)
            , paste0("sigma", 1:Nk)
            , paste0("lambda", 1:Nk)
            , paste0("K", 1:Nv)
            , paste0("Z", 1:Nv))
sims <- mcmc_array(ns, nchains, params)
for (ic in 1:nchains) {
  res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)
  sims[, ic, ] <- cbind(res$mu, res$sigma, res$lambda, res$K, res$Z)
}