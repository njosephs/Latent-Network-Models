library(twitteR)
setup_twitter_oauth("kZADq17e01p6ILN0u3vdgU1DH"
                    , "7bKnY3aGD1LZIDWbS6tVzru8zFZNV492oN9I51g7b6HctxC8Kb"
                    , "266449908-riSJfN9LFhHcRlm2NNkaA2EA30JGHEZf5rBnJPZR"
                    , "BZVZx7nxSAYbsYkpCnPQLpVVlQyFMFIaZ2T49vSbPFasE")

# BU ID
BU <- getUser("BU_Tweets")
# Get BU followers
followers.BU <- BU$getFollowers(n=NULL)
# Get (public) followers' followers
followers_public.BU <- sapply(followers.BU, FUN = function(x) {
  follower <- getUser(x)
  if (follower$protected == TRUE) { # check if public
    next
  } else {
    follower$getFollowers(n = NULL) # if public, get followers
  }})