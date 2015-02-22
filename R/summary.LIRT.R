summary.LIRT <- function(LIRT.out){
	
	print("Summary function not yet written.")
	
	# ntopics <- dim(LDAIRT.out$word.chain)[1]
	# nterms <- dim(LDAIRT.out$word.chain)[2]
	# nvotes <- dim(LDAIRT.out$discrim.chain)[1]
	# nvoters <- dim(LDAIRT.out$ideal.chain)[2]
	
	# # Calculate relative evidence provided by each word with respect to topic
	# if (dim(LDAIRT.out$word.chain)[3] > 1) phi.est <- apply(LDAIRT.out$word.chain,c(1,2),mean) else phi.est <- as.matrix(LDAIRT.out$word.chain[,,1])
	# wordtopicpropensity <- phi.est/t(matrix(rep(colSums(as.matrix(phi.est)),ntopics), nterms ,ntopics))
	# wordtopicevidence <- phi.est*wordtopicpropensity
	
	# # Calculate voter orderings in each topic
	# if (dim(LDAIRT.out$ideal.chain)[3] > 1) theta.est <- apply(LDAIRT.out$ideal.chain,c(1,2),mean) else theta.est <- as.matrix(LDAIRT.out$ideal.chain[,,1])
	# rank.est <- t(apply(theta.est,1,rank))
	
	# # Find votes most heavily in each topic
	# if (dim(LDAIRT.out$topic.chain)[3] > 1) lambda.est <- apply(LDAIRT.out$topic.chain,c(1,2),mean) else lambda.est <- as.matrix(LDAIRT.out$topic.chain[,,1])

	# TopWords <- matrix(NA,20,ntopics)
	# TopVotes <- matrix(NA,20,ntopics)
	
	# for (k in 1:ntopics){
		# print(paste("Topic ",k,":",sep=""))
		# TopWords[,k] <- sort(wordtopicevidence[k,],decreasing=TRUE,index.return=TRUE)$ix[1:20]
		# print(paste("Distinctive Words: ",paste(names(sort(wordtopicevidence[k,],decreasing=TRUE))[1:10],collapse=", "),sep=""))
		# TopVotes[,k] <- sort(lambda.est[k,],decreasing=TRUE,index.return=TRUE)$ix[1:20]
		# print(paste("Most Topical Votes: ",paste(paste(colnames(lambda.est)[sort(lambda.est[k,],decreasing=TRUE,index.return=TRUE)$ix],paste("(",round(lambda.est[k,],2),")",sep="")[sort(lambda.est[k,],decreasing=TRUE,index.return=TRUE)$ix])[1:10],collapse="; "),sep=""))
		# print(paste("Voter Ordering: ",paste(names(sort(rank.est[k,],decreasing=FALSE)),collapse=", "),sep=""))
		# # print(paste("Locations: ",paste(sort(theta.est[k,],decreasing=FALSE),collapse=", "),sep=""))
		# print("")
	# }



	# return(list(theta.est=theta.est,lambda.est=lambda.est,phi.est=phi.est,TopWordsByTopic=TopWords,TopVotesByTopic=TopVotes,Ntopics=ntopics,Nterms=nterms,Nvotes=nvotes,Nvoters=nvoters))

}
