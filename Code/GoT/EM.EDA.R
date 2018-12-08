#----------------------------------------
#
#             EM EDA
#
#----------------------------------------

#load up files + libraries
library(ggplot2)
library(gplots)
load("./Data/A.Rdata")
source("./Code/GoT/LNM.EM.R")

#run EM
em <- LNM.EM(A3)

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

D <- matrix(NA, nrow = nrow(A3), ncol = ncol(A3))
D[upper.tri(D)] <- em$d
D[lower.tri(D)] <- t(D)[lower.tri(D)] 
diag(D) <- 0
rownames(D) <- colnames(D) <- rownames(A3)

G <- graph_from_adjacency_matrix(D, 
                                 weighted = TRUE, 
                                 mode = "undirected", 
                                 add.rownames = TRUE)

V(G)$label.cex <- degree(G) / max(degree(G))
plot(G 
     , vertex.size = strength(G) / 10
     , edge.width = log(E(G)$weight)
     #, layout = layout.circle(G)
     , layout = layout_with_dh(G)
     , color = "grey86"
     , vertex.color = "lightgreen"
     , curved = 200)

P <- matrix(NA, nrow = nrow(A3), ncol = ncol(A3))
P[upper.tri(P)] <- em$p
P[lower.tri(P)] <- t(P)[lower.tri(P)] 
diag(P) <- 0
rownames(P) <- colnames(P) <- rownames(A3)

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

