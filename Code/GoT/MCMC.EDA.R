library(rgl)
library(igraph)

load("~/Documents/GitHub/Latent-Network-Models/Data/A.Rdata")
source("~/Documents/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

mcmc_array <- function (ns, nchains, params) {
  nparams <- length(params)
  array(dim = c(ns, nchains, nparams),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain", 1:nchains),
                        parameters = params))
}

ns <- 50000
burn <- 10000
Nv <- nrow(A)
Nk <- 3
d <- 2


set.seed(589)
res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)
samp <- seq(burn, ns, 10) #thinning

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
clusters <- apply(res$K[, samp], 1, Mode)

#######################
# d = 1
#######################

Z.map <- apply(res$Z[, 1, samp], 1, Mode)
mu.map <- apply(res$mu[, 1, samp], 1, mean)
sigma.map <- apply(res$sigma[, samp], 1, mean)

pdf("./Desktop/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/latent_embedding.pdf")
hist(Z.map, freq = FALSE)
for (i in 1:Nk) {
  curve(expr = dnorm(x, mean = mu.map[i], sd = sqrt(sigma.map[i]))
        , col = i + 1, from = -3, to = 3, add = TRUE)  
}
dev.off()

#pdf("./Desktop/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/K_clusters.pdf")
G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G 
     , vertex.color = adjustcolor(clusters + 1, alpha = .5)
     , vertex.label.color = "black"
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
dev.off()

#######################
# d = 2
#######################

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
mu.map <- sapply(1:Nk, FUN = function(i) apply(res$mu[i, , samp], 1, mean))
sigma.map <- apply(res$sigma[, samp], 1, mean)

#pdf("~/Documents/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/latent_embedding.pdf")
plot(Z.map[2, ] ~ Z.map[1, ], col = clusters + 1)
points(t(mu.map), pch = 3, col = 1:Nk + 1)
symbols(t(mu.map), circles = 3*sqrt(sigma.map), lty = 2, fg = 1:Nk + 1, add = TRUE)
dev.off()

#pdf("~/Documents/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/K_clusters.pdf")
G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G 
     , vertex.color = adjustcolor(clusters + 1, alpha = .5)
     , vertex.label.color = "black"
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
dev.off()

# cluster by distance between nodes
cluster.conf <- sapply(1:Nv, FUN = function(i) {
  sum(res$K[i, samp] == clusters[i]) / length(samp)
})

plot(G
     , vertex.color = adjustcolor(clusters + 1, alpha = .5)
     , vertex.label.color = "black"
     , vertex.size = 15 * cluster.conf / min(cluster.conf)
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)


plot(Z.map[2, ] ~ Z.map[1, ], col = clusters + 1, cex = 4*cluster.conf)
points(t(mu.map), pch = 3, col = 1:Nk + 1)
#symbols(t(mu.map), circles = 3*sqrt(sigma.map), lty = 2, fg = 1:Nk + 1, add = TRUE)

#par(mfrow = c(2, 3))  ### here we test the convergence of K for each node
#for (v in 1:Nv) {
#  test <- matrix(nrow = length(samp), ncol = 1)
#  for (i in 1:length(samp)) {
#    test[i, 1] <- sum(res$K[v, samp[1:i]] == clusters[v]) / i
#  }
#  plot(test, type = "l", main = colnames(A)[v])
#  abline(h = 1/3)
#}

#######################
# d = 3
#######################

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
clusters <- kmeans(t(Z.map), Nk)$cluster
plot3d(t(Z.map), col = clusters)

G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G, vertex.color = clusters)