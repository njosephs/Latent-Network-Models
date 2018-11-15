

#INPUT: List of BU public followers, followers_public.BU

adj.list <- vector("list", length = length(followers_public.BU))

for f in followers_public.BU do: 
	 F <- f.getFollowers(n = NULL) #get followers
	 F_in_BU <- intersect(followers_public, F) #subset with public BU followers
	 adj.list[[counter]] <- F_in_BU #add neighborhood to adj_list
	 counter <- counter + 1 #update counter

#OUTPUT: adj.list - the adjcancey list representation of the network
