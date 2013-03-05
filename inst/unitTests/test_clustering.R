

test_dj.jarvisPatrick <- function() {

	numNbrs=10 #5
	minNbrs=2 # 3
	data(sdfsample)
	aps=sdf2ap(sdfsample)
	fp=desc2fp(aps)

#	clustering = jarvisPatrick(aps,j=numNbrs,k=minNbrs,cutoff=0.5,mode="a1a2b",linkage="average")
#	print(clustering)
#	print(table(clustering))
#	sizes=table(clustering)
#	print(sizes[sizes!=1])
#	return()



	for(data in list(aps,fp))
		for(cutoff in list(NULL,0.5))
		{
			nnm = nearestNeighbors(data,numNbrs=numNbrs,cutoff=cutoff)
			for(mode in c("a1a2b","a1b","b"))
				for(linkage in c("single","average","complete")){
					print(paste("cutoff:",cutoff,"mode:",mode,"linkage:",linkage))
					#clustering = jarvisPatrick(data,j=numNbrs,k=minNbrs,cutoff=cutoff,mode=mode,linkage=linkage)
					clustering = jarvisPatrick(nnm,k=minNbrs,mode=mode,linkage=linkage)
					#print(clustering)
					#print(table(clustering))
					sizes=table(clustering)
					print(sizes[sizes!=1])
				}
		}
			
	#clustering = jarvisPatrick(fp,j=numNbrs,k=minNbrs)
	#print(clustering)
}
