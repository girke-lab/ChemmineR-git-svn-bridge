
testDir = system.file("unitTests",package="ChemmineR")
genRandFp = function(n,type=NULL) {
	if(is.null(type))
		new("FP",fp=sapply(runif(n),function(x) if(x>0.5) 1 else 0))
	else
		new("FP",type=type,fp=sapply(runif(n),function(x) if(x>0.5) 1 else 0))
}

test.fp <- function(){
										 

	fp1 = genRandFp(16)
	fp2 = genRandFp(24)
	fp3 = genRandFp(16)

	checkTrue(inherits(fp1,"FP"))

	x=fold(fp1)
	checkEquals(numBits(x),8)
	checkEquals(foldCount(x),1)

	x=fold(fp2,bits=6,count=34) #count should be ignored here
	checkEquals(numBits(x),6)
	checkEquals(foldCount(x),2)

	x=fold(fp3,count=2)
	checkEquals(numBits(x),4)
	checkEquals(foldCount(x),2)

	fp4 = genRandFp(16)
	x=fold(fp4,16)
	checkEquals(numBits(x),1)
	checkEquals(foldCount(x),4)

	x1=genRandFp(16,"test")
	x2=genRandFp(16,"test")
	x3=genRandFp(18,"test")
	x4=genRandFp(16,"test")

	checkException(c(x1,x3))   #different number of bits
	checkException(c(x1,fp1))  #different label
	checkException(c(x1,"hi")) #mixed types

	fpset1 = c(x1,x2,x4)
	checkEquals(length(fpset1),3)

}

test.fpset <- function(){

	data = matrix(replicate(1280,if(runif(1)>0.5)1 else 0),10,128)
    options(bigmemory.typecast.warning=FALSE)
	fpset = new("FPset",fpma=as.big.matrix(data, type="char"))

	checkTrue(inherits(fpset,"FPset"))

	x=fold(fpset)
	checkEquals(numBits(x),64)
	checkEquals(foldCount(x),1)

	x=fold(fpset,bits=16)
	checkEquals(numBits(x),16)
	checkEquals(foldCount(x),3)
	
	x=fold(fpset,count=2)
	checkEquals(numBits(x),32)
	checkEquals(foldCount(x),2)

	x=fold(fpset,16)
	checkEquals(numBits(x),1)
	checkEquals(foldCount(x),7)

	
	fpset2 = new("FPset",fpma=as.big.matrix(data, type="char"),type="randFP")
	checkEquals(fptype(fpset2),"randFP")
	checkException(c(fpset,fpset2))

	data = matrix(replicate(640,if(runif(1)>0.5)1 else 0),10,64)
	fpset3 = new("FPset",fpma=as.big.matrix(data, type="char"))
	checkException(c(fpset,fpset3))

	fpset4 = new("FPset",fpma=as.big.matrix(data, type="char"),type="randFP")
	catfp1=c(fpset2,fpset2,fpset2)
	checkEquals(length(catfp1),30)

}

test.import <- function(){

	#test single molecule sdf files
	singleMolFile = file.path(testDir,"singleMolecule.sdf")
	singleMols = read.SDFset(singleMolFile)
	checkEquals(length(singleMols),2)

	#test mol file import
	checkEquals(length(read.SDFset(file.path(testDir,"sample.mol"))),1)

	#test sdf file import
	checkEquals(length(read.SDFset(file.path(testDir,"sample.sdf"))),1)

}
