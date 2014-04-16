


test.formatConversions <- function() {
	message("test.formatConversions")

	data(sdfsample)

	smiles = sdf2smiles(sdfsample[1:10])
	sdfs = smiles2sdf(smiles)
	smiles2 = sdf2smiles(sdfs)

	checkEquals(smiles,smiles2)

}

test.genAPDescriptors <- function(){

	data(sdfsample)
	
	for(i in 1:100){
		sdf = sdfsample[[i]]
		desc = genAPDescriptors(sdf)
		#print(desc);

		oldDesc=ChemmineR:::.gen_atom_pair(ChemmineR:::SDF2apcmp(sdf))
                 
		#print(oldDesc);
		checkTrue(all(desc == oldDesc))
		#print(all(desc == oldDesc))
	}

}
test.propOB <- function() {
	data(sdfsample)
	p = propOB(sdfsample[1:5])
	#print(p)
	#checkEquals(ncol(p),15)
	checkEquals(nrow(p),5)
   checkEquals(p$MW[2],MW(sdfsample[2])[[1]])

}
test.fingerprintOB <- function(){
	if(require(ChemmineOB)){
		data(sdfsample)
		fp = fingerprintOB(sdfsample[1:5],"FP2")
		print(fp)
		checkEquals(fptype(fp),"FP2")
	}
}
test.obmolRefs <- function() {
	data(sdfsample)
	if(require(ChemmineOB)){
		obmolRef = obmol(sdfsample[[1]])
		checkEquals(class(obmolRef),"_p_OpenBabel__OBMol")

		obmolRefs = obmol(sdfsample)
		checkEquals(class(obmolRefs),"list")
		checkEquals(class(obmolRefs[[2]]),"_p_OpenBabel__OBMol")
		checkEquals(length(sdfsample),length(obmolRefs))
	}else
		checkException(obmol(sdfsample[[1]]))

}
test.smartsSearchOB <- function(){
	data(sdfsample)
	rotableBonds = smartsSearchOB(sdfsample[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)
	print("rotable bonds: ")
	print(rotableBonds)
	print(sdfid(sdfsample[1:5]))
	checkEquals(as.vector(rotableBonds[1:5]),c(24,20,14,30,10))

}
test.fpSim <- function(){
	data(apset)
	fpset = desc2fp(apset)
	dists = fpSim(fpset[[1]],fpset,top=6)
	checkEqualsNumeric(dists, 
							 c(1.0000000,0.4719101,0.4288499,0.4275229,0.4247423,0.4187380),
							 tolerance = 0.0001)

	for(m in c("tanimoto","euclidean","tversky","dice")){
		sim = ChemmineR:::fpSimOrig(fpset[[1]],fpset,
					method=m,cutoff=0.4,top=6)
		simFast= fpSim(fpset[[1]],fpset,
					method=m,cutoff=0.4,top=6)
		#message("method: ",m)
		#print(sim)
		#print(simFast)
		checkEqualsNumeric(sim,simFast,tolerance=0.00001)
	}
}
