#----------------------------------------
#
#             EM EDA
#
#----------------------------------------

A <- load("./Data/A.Rdata")



em <- LNM.EM(A)

knitr::kable(data.frame(Alpha = em$alpha, Beta = em$beta))

df <- data.frame(
  x = seq(0, 1, length.out = 1000), 
  y = dgamma(seq(0, 1, length.out = 1000), shape=em$alpha, scale=em$beta))



# no vertex level specific information 
# network level inference on distribution of probs of edges 
# 

