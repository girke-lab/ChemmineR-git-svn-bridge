#!/usr/bin/env Rscript

library(ChemmineR)


#data(sdfsample)
#sdfset = sdfsample

sdfFile = commandArgs(trailingOnly=TRUE)

fpsetDataName = paste(sdfFile,"-fpset",".Rdata",sep="")
sdfsetDataName = paste(sdfFile,".Rdata",sep="")

if(!file.exists(fpsetDataName)){
	if(!file.exists(sdfsetDataName)){
		sdfset = read.SDFset(sdfFile)
		message("done loading")
		sdfset = sdfset[validSDF(sdfset)]
		#sdfset = smiles2sdf(sdf2smiles(sdfset))
		save(sdfset,file=sdfsetDataName)
	}
	else
		load(sdfsetDataName)

	fpset = desc2fp(sdf2ap(sdfset))
	message("desc time: ",ChemmineR:::times$descT)
	message("fac time: ",ChemmineR:::times$facT)
	message("vec time: ",ChemmineR:::times$vecT)

	save(fpset,file=fpsetDataName)

}else
	load(fpsetDataName)



eval(fpset) #force evaluation
fpset=c(fpset,fpset,fpset)
message("starting test of ",length(fpset), " FPs")
print(system.time(ChemmineR:::fpSimOrig(fpset[[1]],fpset)))
#print(system.time(fpSim(fpset[[1]],fpset)))

