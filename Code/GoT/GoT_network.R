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
set.seed(102)
library(igraph)
G <- graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected", add.rownames = TRUE)
V(G)$label.cex <- degree(G) /(max(degree(G)))

pdf("./Final Report/report_figures/initial_net.pdf")
plot(G
     , edge.width = log(E(G)$weight)/mean(log(E(G)$weight))
     , vertex.size = degree(G)
     , vertex.label.color = "black"
     , vertex.color = adjustcolor("lightgreen", alpha.f = .75)
     , layout = layout_with_dh(G)
     , color = adjustcolor("lightgreen", alpha.f = .25)
     , curved = 200)
dev.off()

# Only keep important characters
A2 <- A[which(rowSums(A) > 100), which(rowSums(A) > 100)]
G2 <- graph_from_adjacency_matrix(A2, weighted = TRUE, mode = "undirected", add.rownames = TRUE)
V(G2)$label.cex <- degree(G2) /(max(degree(G2)))

pdf("./Final Report/report_figures/filtered_net.pdf")
plot(G2
     , edge.width = log(E(G2)$weight)/mean(log(E(G2)$weight))
     , vertex.size = 2*degree(G2)
     , vertex.label.color = "black"
     , vertex.color = adjustcolor("lightgreen", alpha.f = .75)
     , layout = layout_with_dh(G2)
     , color = adjustcolor("lightgreen", alpha.f = .25)
     , curved = 200)
dev.off()
#final unweighted/weighted adjacency matrix
A <- ifelse(A2 != 0, 1, 0)
W <- A2

#save Adjacency matrices 
save(A, file = "./Data/A.Rdata")
save(W, file = "./Data/W.Rdata")

