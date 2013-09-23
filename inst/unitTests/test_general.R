


test.formatConversions <- function() {

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
	if(require(ChemmineOB)){
		data(sdfsample)
		p = propOB(sdfsample[1:5])
		print(p)
		checkEquals(ncol(p),14)
		checkEquals(nrow(p),5)
		checkEquals(p$MW[2],MW(sdfsample[2])[[1]])
	}

}
