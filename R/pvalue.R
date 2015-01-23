
genParameters = function(fpset,distance = fpSim,... ) {
	if( ! inherits(fpset,"FPset"))
		stop("fpset must be an instance of FPset")
	if(length(fpset) == 0)
		stop("fpset cannot be empty")

	totalBitCount = numBits(fpset[1])
	N=length(fpset)

	size=totalBitCount+2
	parameters = data.frame(count=rep(0,size),avg=rep(0,size),varience=rep(0,size),
									alpha=rep(0,size),beta=rep(0,size))

	setBitCounts = apply(as.matrix(fpset),c(1),function(fp) sum(fp))
	stats = data.frame(count=rep(0,size),sums=rep(0,size),squares=rep(0,size))

	for(i in seq(along=fpset)){
		distances = distance(fpset[[i]],fpset,...)
		numBitsSet= setBitCounts[i]+1
		stats$count[numBitsSet]= stats$count[numBitsSet]+ length(distances)
		stats$sums[numBitsSet]= stats$sums[numBitsSet]+ sum(distances)
		stats$squares[numBitsSet]= stats$squares[numBitsSet]+ sum(distances^2)

		#for quries with bit counts with no stats

		#message(numBitsSet," ",stats$count[numBitsSet]," ",length(distances))
		stats$count[size]= stats$count[size]+ stats$count[numBitsSet]
		stats$sums[size]= stats$sums[size]+ stats$sums[numBitsSet]
		stats$squares[size]= stats$squares[size]+ stats$squares[numBitsSet]
	}
	message("stats: ")
	print(stats[40,])
	print(stats[367,])

	for(i in 1:(totalBitCount+2)){
		avg =stats$sums[i] / stats$count[i]
		varience = (stats$squares[i] - (stats$sums[i]^2/stats$count[i]))/stats$count[i]

		parameters$count[i] = stats$count[i]
		parameters$avg[i] = avg
		parameters$varience[i]=varience
		parameters$alpha[i] = avg*((avg*(1-avg))/varience - 1)
		parameters$beta[i] = (1-avg)*((avg*(1-avg))/varience -1)
	}

	parameters
}

scoredDistance = function(query,db,parameters,distance=fpSim,...){

	distances= distance(query,db,...)

	numBitsSet = sum(as.numeric(query))+1

	N = parameters$count[numBitsSet]
	if(N==0){ # no stats collected for this number of bits so use global values
		warning("no parameters avaliable for fingerprints with ",numBitsSet-1," bits set, using global parameters")
		numBitsSet=nrow(parameters) #global stats are last element of parameters
		N = parameters$count[numBitsSet]
	}

	message("using stats for ",numBitsSet-1," bits")
	print(parameters[numBitsSet,])
	
   avg = parameters$avg[numBitsSet]
	varience = parameters$varience[numBitsSet]
	alpha = parameters$alpha[numBitsSet]
	beta = parameters$beta[numBitsSet]

	Reduce(rbind,Map(function(i){
			 zscore = (distances[i] - avg) / sqrt(varience)
			 evalue = N*(1-pbeta(distances[i],alpha,beta))
			 pvalue = 1-exp(-evalue)
			 data.frame(distance=distances[i],zscore=zscore,evalue=evalue,pvalue=pvalue)
		},seq(along=distances)))

}
