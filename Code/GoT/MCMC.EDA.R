library(igraph)
# library(mixtools) # for 2D contour plots
# library(rgl)      # for 3D plots

# load("~/Documents/GitHub/Latent-Network-Models/Data/A.Rdata")
# source("~/Documents/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")
load("./Desktop/GitHub/Latent-Network-Models/Data/A.Rdata")
source("./Desktop/GitHub/Latent-Network-Models/Code/GoT/LNM.MCMC.R")

ns <- 10000
burn <- 4000
Nv <- nrow(A)
Nk <- 4
d <- 4

set.seed(589)
res <- LNM.MCMC(A, Nk = Nk, d = d, ns = ns)
samp <- seq(burn, ns, 10) # thinning

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
clusters <- apply(res$K[, samp], 1, Mode)
cluster.conf <- sapply(1:Nv, FUN = function(i) {
  sum(res$K[i, samp] == clusters[i]) / length(samp)
})

par(mfrow = c(2, 3))
for (v in 1:Nv) {
  test <- matrix(nrow = length(samp), ncol = 1)
  for (i in 1:length(samp)) {
    test[i, 1] <- sum(res$K[v, samp[1:i]] == clusters[v]) / i
  }
  plot(test, type = "l", main = colnames(A)[v])
  abline(h = 1/Nk)
}  
dev.off()

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
     , vertex.size = 15 * cluster.conf / min(cluster.conf)
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
plot(Z.map[2, ] ~ Z.map[1, ], col = clusters + 1, cex = Nk*cluster.conf)
points(t(mu.map), pch = 3, col = 1:Nk + 1)
# plot(Z.map[2, ] ~ Z.map[1, ], col = clusters + 1, cex = 4*cluster.conf, xlim = c(-1.5, 1.5), ylim = c(-2, 2))
# points(t(mu.map), pch = 3, col = 1:Nk + 1)
# for (i in 1:Nk) {
#   # proper CI
#   ellipse(t(mu.map)[i, ], sigma.map[i]*diag(d), col = i + 1, alpha = .33, draw = TRUE)
# }
# # incorrect CI
# symbols(t(mu.map), circles = 3*sqrt(sigma.map), lty = 2, fg = 1:Nk + 1, add = TRUE)
dev.off()

#pdf("~/Documents/GitHub/Latent-Network-Models/Final Report/report_figures/MCMC/K_clusters.pdf")
G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G 
     , vertex.color = adjustcolor(clusters + 1, alpha = .5)
     , vertex.label.color = "black"
     , vertex.size = 15 * cluster.conf / min(cluster.conf)
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
dev.off()

#######################
# d = 3+
#######################

Z.map <- sapply(1:Nv, FUN = function(i) apply(res$Z[i, , samp], 1, mean))
clusters <- kmeans(t(Z.map), Nk)$cluster
# plot3d(t(Z.map), col = clusters)

G <- graph_from_adjacency_matrix(A
                                 , mode = "undirected"
                                 , add.rownames = TRUE)
plot(G 
     , vertex.color = adjustcolor(clusters + 1, alpha = .5)
     , vertex.label.color = "black"
     , vertex.size = 15 * cluster.conf / min(cluster.conf)
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
