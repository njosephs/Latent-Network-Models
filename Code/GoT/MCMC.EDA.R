load("./Desktop/GitHub/Latent-Network-Models/Data/A.Rdata")
source("./Desktop/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

ns <- 50000
burn <- 20000
Nv <- nrow(A)
Nk <- 4
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
Z.map <- apply(res$Z[, 1, (burn+1):ns], 1, mean)
hist(Z.map)
plot(density(Z.map))
plot(sort(Z.map))
clusters <- kmeans(Z.map, Nk)$cluster
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]
colnames(A)[clusters == 4]
colnames(A)[clusters == 5]

##### d = 2
# plot(res$Z[4, 1, ], type = "l")

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , (burn+1):ns], 1, mean))
plot(Z.map[1, ] ~ Z.map[2, ])

clusters <- kmeans(t(Z.map), 2)$cluster
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]
colnames(A)[clusters == 4]

# ns = 50000, burn = 20000, 2 clusters
# Starks vs

mu.map <- sapply(1:Nk, FUN = function(i) apply(res$mu[i, , (burn+1):ns], 1, mean))
points(t(mu.map), pch = 3)

##### Nk = 2
K.probs <- apply(res$K[, ] - 1, 1, mean)
plot(sort(K.probs))

Z.conf <- ifelse(K.probs < .495, 1, ifelse(K.probs > .505, 2, "unsure"))
colnames(A)[Z.conf == 1]
colnames(A)[Z.conf == 2]
colnames(A)[Z.conf == "unsure"]
