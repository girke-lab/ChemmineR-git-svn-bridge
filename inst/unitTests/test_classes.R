
genRandFp = function(n) new("FP",fp=sapply(runif(n),function(x) if(x>0.5) 1 else 0))
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

}

test.fpset <- function(){

	data = matrix(replicate(1280,if(runif(1)>0.5)1 else 0),10,128)
	fpset = new("FPset",fpma=data)

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

	
	fpset2 = new("FPset",fpma=data,type="randFP")
	checkEquals(fptype(fpset2),"randFP")
	checkException(c(fpset,fpset2))

	data = matrix(replicate(640,if(runif(1)>0.5)1 else 0),10,64)
	fpset2 = new("FPset",fpma=data)
	checkException(c(fpset,fpset2))
}
