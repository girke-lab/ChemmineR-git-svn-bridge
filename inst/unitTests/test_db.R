
test_aaa.clean <- function(){

	unlink(c("test.db","test.sdf"))
	
}

test_aa.initDb<-function(){
	conn = initDb("test.db")
	checkTrue(file.exists("test.db"))
}

test_ba.loadSdf<-function(){
	data(sdfsample)
	tempFile=tempfile()
	write.SDF(sdfsample,tempFile)
	conn = initDb("test.db")
	loadSdf(conn,tempFile,function(sdfset)cbind(MW=MW(sdfset)) )

	compoundCount = dbGetQuery(conn,"SELECT count(*) FROM compounds")[1][[1]]
	checkEquals(compoundCount ,length(cid(sdfsample)))
	featureCount= dbGetQuery(conn,"SELECT count(*) FROM feature_MW")[1][[1]]
	checkEquals(featureCount ,length(cid(sdfsample)))
}

test_ca.findCompounds<-function(){

	conn = initDb("test.db")

	indexes = findCompounds(conn,"MW","MW < 400")
	print(paste("found",length(indexes)," compounds"))

}

test_da.getCompounds<-function(){

	conn = initDb("test.db")

	indexes = findCompounds(conn,"MW","MW < 400")

	sdfset = getCompounds(conn,indexes)
	checkEquals(length(cid(sdfset)),70)
	
	getCompounds(conn,indexes,file="test.sdf")
	checkTrue(file.exists("test.sdf"))
	sdfFromFile = read.SDFset("test.sdf")
	checkEquals(length(cid(sdfFromFile)),70)

}
