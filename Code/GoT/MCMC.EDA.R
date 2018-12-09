library(rgl)
library(igraph)

# load("~/Documents/GitHub/Latent-Network-Models/Data/A.Rdata")
# source("~/Documents/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

load("./Desktop/GitHub/Latent-Network-Models/Data/A.Rdata")
source("./Desktop/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

ns <- 25000
burn <- 10000
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
samp <- seq(burn, ns, 10) # thinning

##### d = 2
# plot(res$Z[4, 1, ], type = "l")

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
clusters <- apply(res$K[, samp], 1, Mode)
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]
plot(Z.map[1, ] ~ Z.map[2, ], col = clusters)
mu.map <- sapply(1:Nk, FUN = function(i) apply(res$mu[i, , samp], 1, mean))
points(t(mu.map), pch = 3)
sigma.map <- apply(res$sigma[, samp], 1, mean)
symbols(t(mu.map), circles = sigma.map, lty = 2, add = TRUE)

G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G, vertex.color = clusters)

##### d = 3
Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
clusters <- kmeans(t(Z.map), Nk)$cluster
plot3d(t(Z.map), col = clusters)

G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G, vertex.color = clusters)

##### d = 1
Z.map <- apply(res$Z[, 1, samp], 1, mean)
hist(Z.map)
plot(density(Z.map))
plot(sort(Z.map))
clusters <- kmeans(Z.map, Nk)$cluster
colnames(A)[clusters == 1]
colnames(A)[clusters == 2]
colnames(A)[clusters == 3]
colnames(A)[clusters == 4]
colnames(A)[clusters == 5]

# ns = 50000, burn = 20000, thin by 10, 2 clusters
# Starks vs

##### Nk = 2
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
K.probs <- apply(res$K[, ] - 1, 1, Mode)
plot(sort(K.probs))

Z.conf <- ifelse(K.probs < .495, 1, ifelse(K.probs > .505, 2, "unsure"))
colnames(A)[Z.conf == 1]
colnames(A)[Z.conf == 2]
colnames(A)[Z.conf == "unsure"]