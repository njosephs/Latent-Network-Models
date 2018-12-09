#----------------------------------------
#
#             EM Unweighted EDA
#
#----------------------------------------

#load up files + libraries---------------
library(ggplot2)
library(gplots)
load("./Data/A.Rdata")
source("./Code/GoT/LNM.EM.R")

set.seed(1985)

#run EM
em <- LNM.EM.U(A)

#Plot Networks
D <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
D[upper.tri(D)] <- em$d
D[lower.tri(D)] <- t(D)[lower.tri(D)] 
diag(D) <- 0
rownames(D) <- colnames(D) <- rownames(A)

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2 *max(degree(G)))
layout <- layout_with_dh(G)

pdf("./Final Report/report_figures/graph_dist_unweighted.pdf")
plot(G 
     , vertex.size = degree(G)
     , edge.width = log(E(G)$weight)
     , layout = layout
     , color = "grey86"
     , vertex.color = adjustcolor("green", alpha.f = .75)
     , curved = 200)
dev.off()

P <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
P[upper.tri(P)] <- exp(em$pi)
P[lower.tri(P)] <- t(P)[lower.tri(P)] 
diag(P) <- 0
rownames(P) <- colnames(P) <- rownames(A)

G <- graph_from_adjacency_matrix(P, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2*max(degree(G)))

pdf("./Final Report/report_figures/graph_p_unweighted.pdf")
plot(G 
     , vertex.size = degree(G)
     , edge.width = E(G)$weight
     , layout = layout
     , color = "grey86"
     , vertex.color = adjustcolor("red", alpha.f = .75)
     , curved = 200)
dev.off()

#Plot heatmps
df <- melt(D)

pdf("./Final Report/report_figures/heatmap_dist_unweighted.pdf")
ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "white", mid = "red", midpoint = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "Distance")
dev.off()

df <- melt(P)

pdf("./Final Report/report_figures/heatmap_p_unweighted.pdf")
ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "Probabilities")
dev.off()


#----------------------------------------
#
#             EM Weighted EDA
#
#----------------------------------------

#load data 
load("./Data/W.Rdata")

#run EM
em <- LNM.EM.W(W)

#plot density
plot(density(em$d[em$d < 5]))

D <- matrix(NA, nrow = nrow(W), ncol = ncol(W))
D[upper.tri(D)] <- em$d
D[lower.tri(D)] <- t(D)[lower.tri(D)] 
diag(D) <- 0
rownames(D) <- colnames(D) <- rownames(W)
df <- melt(D)

P <- matrix(NA, nrow = nrow(W), ncol = ncol(W))
P[upper.tri(P)] <- em$p
P[lower.tri(P)] <- t(P)[lower.tri(P)] 
diag(P) <- 0
rownames(P) <- colnames(P) <- rownames(W)
df <- melt(P)
p1 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()


df <- melt(W - P)
p2<- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0) 
p2
grid.arrange(p1,p1, ncol = 1)


#----------------------------------------
#
#             Spectral clustering
#
#----------------------------------------

library(rgl)
spectral_clust <- function(S, d, k = 4){
  sqrtD <- diag(1/sqrt(rowSums(S)))
  Lsym <- diag(rep(1, nrow(S))) - tcrossprod(crossprod(sqrtD, S), sqrtD)
  ES <- eigen(Lsym)
  coord <- (ES$vectors[,order(ES$values)])[,2:(d+1)]
  km <- kmeans(coord, k)
  list(groups = km$cluster, coord = coord)
}

G <- graph_from_adjacency_matrix(P, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)
sc <- spectral_clust(P, d = 3, k = 4)
plot3d(sc$coord, col = sc$groups)
plot(G, vertex.color = sc$groups)

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)
sc <- spectral_clust(D, d = 3, k = 4)
plot3d(sc$coord, col = sc$groups)
plot(G, vertex.color = sc$groups)





