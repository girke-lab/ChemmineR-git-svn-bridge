


test.formatConversions <- function() {

	data(sdfsample)

	smiles = sdf2smiles(sdfsample[1:10])
	sdfs = smiles2sdf(smiles)
	smiles2 = sdf2smiles(sdfs)

	checkEquals(smiles,smiles2)

}

test.genDescriptors <- function(){

	data(sdfsample)
	
	for(i in 1:100){
		sdf = sdfsample[[i]]
		desc = genDescriptors(sdf)
		desc = ChemmineR:::.factor_to_vector(as.factor(desc))
		#print(desc);

		oldDesc = ap(sdf2ap(sdf))
		#print(oldDesc);
		print(all(desc == oldDesc))
	}

}
