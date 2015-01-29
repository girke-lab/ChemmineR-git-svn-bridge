


test.formatConversions <- function() {
	message("test.formatConversions")

	data(sdfsample)

	smiles = sdf2smiles(sdfsample[1:10])
	sdfs = smiles2sdf(smiles)
	smiles2 = sdf2smiles(sdfs)

	checkEquals(smiles,smiles2)

}

test.genAPDescriptors <- function(){

	DEACTIVATED("removed old version of function")

	data(sdfsample)
	
	for(i in 1:100){
		sdf = sdfsample[[i]]
		desc = genAPDescriptors(sdf) 
		#print(head(desc));

		oldDesc=ChemmineR:::.gen_atom_pair(ChemmineR:::SDF2apcmp(sdf))
                 
		#print(oldDesc);
		compResult = desc==oldDesc
		if(!all(compResult)){
			message("descriptor mismatch")
			print(oldDesc[!compResult])
			print(desc[!compResult])
			message("----------")
			firstFalse = match(FALSE,compResult)
			print(oldDesc[(firstFalse-5):(firstFalse+5)])
			print(desc[(firstFalse-5):(firstFalse+5)])

		}
		checkTrue(all(compResult))
		#checkEqualsNumeric(desc,oldDesc)
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
		checkEquals(fptype(fp),"FP2")
		fpSingle = fingerprintOB(sdfsample[1],"FP2")
		checkEquals(as.character(class(fpSingle)),"FPset")
		checkEqualsNumeric(as.matrix(fpSingle[1]), as.matrix(fp[1]))
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
test.exactMassOB <- function(){
	data(sdfsample)
	mass = exactMassOB(sdfsample[1:5])
	checkEqualsNumeric(mass,c(456.2009,357.1801,
									  370.1100,461.1733,
									  318.1943),tolerance=0.00001)
}
test.3dCoords <-function(){
	data(sdfsample)
	sdf3d = generate3DCoords(sdfsample[1])

	checkTrue(!any(atomblock(sdf3d)[[1]][,3]==0))
	
}
test.canonicalize <- function(){
	data(sdfsample)
	cansdf = canonicalize(sdfsample[1])

	bb=bondblock(cansdf)[[1]]

	checkEqualsNumeric(bb[1,1:3],c(2,3,1))
	checkEqualsNumeric(bb[2,1:3],c(2,4,1))

}
test.parseV3000 <- function() {

	DEACTIVATED("requires local files")
	sdfset2 = read.SDFset("~/runs/v3000/DrugLike-0_2-3K3K_1.v2k.sdf")
	sdfset3 = read.SDFset("~/runs/v3000/DrugLike-0_2-3K3K_1.sdf")  

	compareSdfVersions = function(v2k,v3k){
		checkEquals(sdfid(v2k),sdfid(v3k))

		#message("v2k: ",nrow(atomblock(v2k)),"x",ncol(atomblock(v2k)))
		#message("v3k: ",nrow(atomblock(v3k)),"x",ncol(atomblock(v3k)))
		#print(head(atomblock(v2k)))
		#print(head(atomblock(v3k)))
		toCompare = c(1:5,7:10) #exclude colum 6
		checkTrue( all(atomblock(v2k)[,toCompare] ==
							atomblock(v3k)[,toCompare]))

#		cmp = bondblock(v2k)[,1:3] == bondblock(v3k)
#		if(! all(cmp)){
#			mismatched = which(cmp==FALSE)
#			print(cmp)
#			print("mismatched: ")
#			print(mismatched)
#			print("data:")
#			print(bondblock(v2k)[mismatched,1:3])
#			print(bondblock(v3k)[mismatched,])
#		}
#		checkTrue( all(bondblock(v2k)[,1:3] == bondblock(v3k)))

		checkTrue( all(datablock(v2k) == datablock(v3k)))
	}
	for(i in seq(along=sdfset2)){
		#if(!(i  %in% c(38,39,89))){ # this differ in acceptable ways
		#	message("testing ",i, " id: ",sdfid(sdfset2[i]))
			compareSdfVersions(sdfset2[[i]],sdfset3[[i]])
		#}
	}


}
