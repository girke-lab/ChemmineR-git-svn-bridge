

test_dj.jarvisPatrick <- function() {

	numNbrs=5
	minNbrs=2
	data(sdfsample)
	aps=sdf2ap(sdfsample)
	fp=desc2fp(aps)

	for(data in list(aps,fp))
		for(cutoff in list(NA,0.7))
			for(mode in c("a1a2b","a1b","b")){
				print(paste("cutoff:",cutoff,"mode:",mode))
				clustering = jarvisPatrick(data,j=numNbrs,k=minNbrs,cutoff=cutoff,mode=mode)
				counts = table(clustering)
				print(counts[counts > 1])
				#print(clustering)
			}
	#clustering = jarvisPatrick(fp,j=numNbrs,k=minNbrs)
	#print(clustering)
}
