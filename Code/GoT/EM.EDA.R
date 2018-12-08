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
em <- LNM.EM(A)

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



# no vertex level specific information 
# network level inference on distribution of probs of edges 
# 

