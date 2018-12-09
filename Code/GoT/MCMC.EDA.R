library(rgl)
library(igraph)

# load("~/Documents/GitHub/Latent-Network-Models/Data/A.Rdata")
# source("~/Documents/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

load("./Desktop/GitHub/Latent-Network-Models/Data/A.Rdata")
source("./Desktop/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

ns <- 10000
burn <- 2000
Nv <- nrow(A)
Nk <- 5
d <- 2

set.seed(589)
res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)
samp <- seq(burn, ns, 10) # thinning

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

pdf("./Desktop/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/K_clusters.pdf")
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

pdf("./Desktop/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/latent_embedding.pdf")
plot(Z.map[2, ] ~ Z.map[1, ], col = clusters + 1)
points(t(mu.map), pch = 3, col = 1:Nk + 1)
symbols(t(mu.map), circles = 3*sqrt(sigma.map), lty = 2, fg = 1:Nk + 1, add = TRUE)
dev.off()

pdf("./Desktop/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/K_clusters.pdf")
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
# d = 3
#######################

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
clusters <- kmeans(t(Z.map), Nk)$cluster
plot3d(t(Z.map), col = clusters)

G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G, vertex.color = clusters)