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
     , vertex.label.color = "black"
     , edge.width = log(E(G)$weight)
     , layout = layout
     , color = "grey86"
     , vertex.color = adjustcolor("lightblue", alpha.f = .75)
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
     , vertex.label.color = "black"
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
  scale_fill_gradient2(low = "blue", high = "white", mid = "lightblue")+
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

#Density plots
pdf("./Final Report/report_figures/density_dist_weighted.pdf")
ggplot(data.frame(dist = em$d[em$d <5]), aes(x=dist)) + 
  geom_density(color="grey", 
               fill="lightblue", 
               alpha = .7)+
  theme_minimal()+
  labs(x = "Latent Distance", 
       y = "Density", 
       title = "Non-Zero Distance Distribution")
dev.off()

#Network Figures
D <- matrix(NA, nrow = nrow(W), ncol = ncol(W))
D[upper.tri(D)] <- em$d
D[lower.tri(D)] <- t(D)[lower.tri(D)] 
diag(D) <- 0
rownames(D) <- colnames(D) <- rownames(W)

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2 *max(degree(G)))

pdf("./Final Report/report_figures/graph_dist_weighted.pdf")
plot(G 
     , vertex.size = degree(G)
     , vertex.label.color = "black"
     , edge.width = log(E(G)$weight)
     , layout = layout
     , color = "grey86"
     , vertex.color = adjustcolor("lightblue", alpha.f = .75)
     , curved = 200)
dev.off()


P <- matrix(NA, nrow = nrow(W), ncol = ncol(W))
P[upper.tri(P)] <- em$p
P[lower.tri(P)] <- t(P)[lower.tri(P)] 
diag(P) <- 0
rownames(P) <- colnames(P) <- rownames(W)

G <- graph_from_adjacency_matrix(P, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2*max(degree(G)))

pdf("./Final Report/report_figures/graph_p_weighted.pdf")
plot(G 
     , vertex.size = degree(G)
     , vertex.label.color = "black"
     , edge.width = log(E(G)$weight)
     , layout = layout
     , color = "grey86"
     , vertex.color = adjustcolor("red", alpha.f = .75)
     , curved = 200)
dev.off()

df <- melt(D)
p1 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "white", mid = "lightblue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "Distance")

pdf("./Final Report/report_figures/heatmap_dist_weighted.pdf")
p1
dev.off()

df <- melt(P)
p2 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "Mean")

pdf("./Final Report/report_figures/heatmap_p_weighted.pdf")
p2
dev.off()

df <- melt(W - P)
p3 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "W - P")

pdf("./Final Report/report_figures/heatmap_p_diff_weighted.pdf")
p3
dev.off()

#combo plots
library(gridExtra)

df <- melt(W)
p4 <- ggplot(data = df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Character Names",
       y = "Character Names", 
       fill = "Obs. Weight")


pdf("./Final Report/report_figures/heatmap_comp_weighted.pdf")
grid.arrange(p4, p2,ncol = 2)
dev.off()

#----------------------------------------
#
#             Spectral clustering
#
#----------------------------------------
#set up plotly
library(plotly)
Sys.setenv("plotly_username"="kkung")
Sys.setenv("plotly_api_key"="VHGITshBjV5oaugEaP49")

spectral_clust <- function(S, d, k = 4){
  sqrtD <- diag(1/sqrt(rowSums(S)))
  Lsym <- diag(rep(1, nrow(S))) - tcrossprod(crossprod(sqrtD, S), sqrtD)
  ES <- eigen(Lsym)
  coord <- (ES$vectors[,order(ES$values)])[,2:(d+1)]
  km <- kmeans(coord, k)
  list(groups = km$cluster, coord = coord)
}

#P Clustering --------------------------
sc <- spectral_clust(P, d = 3, k = 4)
eigenvec <- data.frame(sc$coord)
colnames(eigenvec) <- c("oned", "twod", "threed")
p <- plot_ly(eigenvec, x = ~oned, y = ~twod, z = ~threed, color = as.factor(sc$groups/4)) %>%
  add_markers() %>%
  layout(title = "Plot of Eigenvectors for Pi",
           scene = list(xaxis = list(title = 'x'),
                      yaxis = list(title = 'y'),
                      zaxis = list(title = 'z')),
         showlegend = FALSE)
p
chart_link = api_create(p, filename="three_d_P")
#chart_link

G <- graph_from_adjacency_matrix(P, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2 *max(degree(G)))

pdf("./Final Report/report_figures/graph_p_clust_weighted.pdf")
plot(G 
     , vertex.color = adjustcolor(sc$groups +1, alpha = .5)
     , vertex.label.color = "black"
     , edge.width = log(E(G)$weight)
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
dev.off()

#D Clustering --------------------------
sc <- spectral_clust(D, d = 3, k = 4)
plot3d(sc$coord, col = sc$groups)
plot(G, vertex.color = sc$groups)

eigenvec<-data.frame(sc$coord)
colnames(eigenvec)<-c("oned", "twod", "threed")
p <- plot_ly(eigenvec, x = ~oned, y = ~twod, z = ~threed,
             color = as.factor(sc$groups/4)) %>%
  add_markers() %>%
  layout(title = "Plot of Eigenvectors for Distance",
         scene = list(xaxis = list(title = 'x'),
                      yaxis = list(title = 'y'),
                      zaxis = list(title = 'z')), 
         showlegend = FALSE) 
p
chart_link = api_create(p, filename="three_d_D")
#chart_link

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <-  degree(G)/(2 *max(degree(G)))

pdf("./Final Report/report_figures/graph_dist_clust_weighted.pdf")
plot(G 
     #, vertex.size = degree(G)
     , vertex.color = adjustcolor(sc$groups +1, alpha = .5)
     , vertex.label.color = "black"
     , edge.width = log(E(G)$weight)
     , edge.color = adjustcolor("grey86", alpha = .75)
     , curved = 200)
dev.off()

