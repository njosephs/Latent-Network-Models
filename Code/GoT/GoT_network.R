#----------------------------------------
#
#           Preprocessing
#
#----------------------------------------

GoT <- read.csv("./Desktop/GitHub/Latent-Twitter-Models/Data/stormofswords.csv")
GoT$Source <- as.character(GoT$Source)
GoT$Target <- as.character(GoT$Target)
GoT$Weight <- as.integer(GoT$Weight)

# Create weighted adjacency matrix
n <- length(unique(c(GoT$Source, GoT$Target))) # number of characters
A <- matrix(0, nrow = n, ncol = n)
colnames(A) <- rownames(A) <- unique(c(GoT$Source, GoT$Target))
for (v in 1:nrow(GoT)) {
  A[GoT[v, "Source"], GoT[v, "Target"]] <- GoT[v, "Weight"]
}
A <- A + t(A)


#----------------------------------------
#
#             EDA
#
#----------------------------------------

library(igraph)
G <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected", add.rownames = TRUE)
V(G)$label.cex <- degree(G) / max(degree(G))
plot(G
     , edge.width = E(G)$weight / max(E(G)$weight)
     , vertex.size = degree(G)
     , layout = layout_with_dh(G)
     , color = "grey86"
     , margin = rep(0, 4))

# Only keep important characters
A2 <- A[which(rowSums(A) > 100), which(rowSums(A) > 100)]
G2 <- graph_from_adjacency_matrix(A2, weighted = TRUE, mode = "undirected", add.rownames = TRUE)
V(G2)$label.cex <- degree(G2) / max(degree(G2))
plot(G2
     , edge.width = log(E(G2)$weight)
     , vertex.size = strength(G2) / 10
     , layout = layout_with_dh(G2)
     , color = "grey86"
     , vertex.color = "lightgreen"
     , margin = rep(-.35, 4)
     , curved = 200)

#final unweighed adjacency matrix
A <- ifelse(A2 > 5 , 1, 0)

#save Adjacency matrix 
save(A, file = "./Code/A.Rdata")

