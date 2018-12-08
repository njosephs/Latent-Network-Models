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
Nk <- 3
d <- 2
# nchains <- 1
# params <- c(paste0("mu", 1:Nk)
#             , paste0("sigma", 1:Nk)
#             , paste0("lambda", 1:Nk)
#             , paste0("K", 1:Nv)
#             , paste0("Z", 1:Nv))
# sims <- mcmc_array(ns, nchains, params)
# for (ic in 1:nchains) {
#   res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)
#   sims[, ic, ] <- cbind(res$mu, res$sigma, res$lambda, res$K, res$Z)
# }

set.seed(589)
res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)

##### d = 1
Z.map <- apply(res$Z[, 1, ], 1, mean)
hist(Z.map)
plot(density(Z.map))
plot(sort(Z.map))
clusters <- kmeans(Z.map, Nk)$cluster
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]

##### d = 2
plot(res$Z[4, 1, ], type = "l")

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , ], 1, mean))
plot(Z.map[1, ] ~ Z.map[2, ])

clusters <- kmeans(t(Z.map), Nk)$cluster
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]

mu.map <- sapply(1:Nk, FUN = function(i) apply(res$mu[i, , ], 1, mean))
sigma.map <- apply(res$sigma, 1, mean)
