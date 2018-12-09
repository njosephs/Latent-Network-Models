#----------------------------------------
#
#             EM Unweighted EDA
#
#----------------------------------------


#load up files + libraries
library(ggplot2)
library(gplots)
load("./Data/A.Rdata")
load("./Data/W.Rdata")
source("./Code/GoT/LNM.EM.R")

#run EM
em <- LNM.EM.U(A)

#EM tables 
knitr::kable(data.frame(Alpha = em$alpha, Beta = em$beta))
table(em$pi)
table(em$d)

#plot posterior distribution
df <- data.frame(
  x = rep(seq(0.01, .99, length.out = 1000),2), 
  y = c(dbeta(seq(0.01, .99, length.out = 1000), shape1=em$alpha, shape2=em$beta),
        dbeta(seq(0.01, .99, length.out = 1000), shape= 1, shape2= 1)),
  group = c(rep("Estimated Distribution", 1000), rep("Prior", 1000))
  )

p1 <- ggplot(df, aes(x = x, y = y, col = group))+
        geom_line()+
        labs(ylab = "Cumulative Beta Density", 
            xlab = "x", 
            title = "EM Esimates for p Density")+
      theme_minimal()
p1

D <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
D[upper.tri(D)] <- em$d
D[lower.tri(D)] <- t(D)[lower.tri(D)] 
diag(D) <- 0
rownames(D) <- colnames(D) <- rownames(A)

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)
V(G)$label.cex <-  strength(G) / max(strength(G))
plot(G 
     #, vertex.size = strength(G) 
     , edge.width = log(E(G)$weight)
     #, layout = layout.circle(G)
     , layout = layout_with_dh(G)
     , color = "grey86"
     , vertex.color = "lightgreen"
     , curved = 200)

P <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
P[upper.tri(P)] <- em$p
P[lower.tri(P)] <- t(P)[lower.tri(P)] 
diag(P) <- 0
rownames(P) <- colnames(P) <- rownames(A)

G <- graph_from_adjacency_matrix(P, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  strength(G) / max(strength(G))
plot(G 
     , vertex.size = 2 * strength(G)
     , edge.width = E(G)$weight
     #, layout = layout.circle(G)
     , layout = layout_with_dh(G)
     , color = "grey86"
     , vertex.color = "lightgreen"
     , curved = 200)

# no vertex level specific information 
# network level inference on distribution of probs of edges 
# 

#----------------------------------------
#
#             EM Weighted EDA
#
#----------------------------------------

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
                       midpoint = 0,) 
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





