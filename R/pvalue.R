
#sample fraction of all pairwise distances
genParameters = function(fpset,similarity= fpSim,sampleFraction=1,... ) {
	if( ! inherits(fpset,"FPset"))
		stop("fpset must be an instance of FPset")
	if(length(fpset) == 0)
		stop("fpset cannot be empty")

	totalBitCount = numBits(fpset[1])
	N=length(fpset)
	sampleSize=floor(N*sqrt(sampleFraction))
	#message("sample size: ",sampleSize,"^2")
	K=0.5  # used to shift the mean towards 0 to stabalize varience computation

	size=totalBitCount+2
	parameters = data.frame(count=rep(0,size),avg=rep(0,size),varience=rep(0,size),
									alpha=rep(0,size),beta=rep(0,size))

	stats = data.frame(count=rep(0,size),sums=rep(0,size),squares=rep(0,size))

	querySetIndecies = sample.int(N,sampleSize)
	targetSetIndecies = sample.int(N,sampleSize)
	targetFpSet = fpset[targetSetIndecies]

	for(i in querySetIndecies){
		distances = similarity(fpset[[i]],targetFpSet,sorted=FALSE,...)
		numBitsSet= sum(as.numeric(fpset[[i]])) + 1

		n=length(distances)
		sums=sum(distances-K)
		squares= sum((distances-K)^2)
		#print(c(i,n,sums,squares,as.numeric(numBitsSet)))

		stats$count[numBitsSet]= stats$count[numBitsSet] + n
		stats$sums[numBitsSet]= stats$sums[numBitsSet] + sums
		stats$squares[numBitsSet]= stats$squares[numBitsSet] + squares

		#for quries with bit counts with no stats

		#message(numBitsSet," ",stats$count[numBitsSet]," ",length(distances))
		stats$count[size]= stats$count[size] + n
		stats$sums[size]= stats$sums[size] + sums
		stats$squares[size]= stats$squares[size] + squares
	}
	#message("stats: ")
	#print(stats[41,])
	#print(stats[368,])

	for(i in 1:(totalBitCount+2)){
		avg =stats$sums[i] / stats$count[i] + K
		varience = (stats$squares[i] - (stats$sums[i]^2/stats$count[i]))/stats$count[i]

		parameters$count[i] = stats$count[i]
		parameters$avg[i] = avg
		parameters$varience[i]=varience
		parameters$alpha[i] = avg*((avg*(1-avg))/varience - 1)
		parameters$beta[i] = (1-avg)*((avg*(1-avg))/varience -1)
	}

	parameters
}

#scoredDistance = function(query,db,parameters,distance=fpSim,...){
#
#	distances= distance(query,db,...)
#
#	numBitsSet = sum(as.numeric(query))+1
#
#	N = parameters$count[numBitsSet]
#	if(N==0){ # no stats collected for this number of bits so use global values
#		warning("no parameters avaliable for fingerprints with ",numBitsSet-1," bits set, using global parameters")
#		numBitsSet=nrow(parameters) #global stats are last element of parameters
#		N = parameters$count[numBitsSet]
#	}
#
#	message("using stats for ",numBitsSet-1," bits")
#	print(parameters[numBitsSet,])
#	
#   avg = parameters$avg[numBitsSet]
#	varience = parameters$varience[numBitsSet]
#	alpha = parameters$alpha[numBitsSet]
#	beta = parameters$beta[numBitsSet]
#
#	Reduce(rbind,Map(function(i){
#			 zscore = (distances[i] - avg) / sqrt(varience)
#			 evalue = N*(1-pbeta(distances[i],alpha,beta))
#			 pvalue = 1-exp(-evalue)
#			 data.frame(distance=distances[i],zscore=zscore,evalue=evalue,pvalue=pvalue)
#		},seq(along=distances)))
#
#}
