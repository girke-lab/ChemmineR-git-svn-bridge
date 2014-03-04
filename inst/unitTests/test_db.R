
test_aaa.clean <- function(){

	unlink(c("test.db","test2.db","test1.db","test.sdf","test_desc.db",
				"part1.db","part2.db","dup_test.db"))
	
}

test_aa.initDb<-function(){
	#DEACTIVATED("temp")
	conn = initDb("test.db")
	checkTrue(file.exists("test.db"))
	checkException(initDb(c(1,2)))
	checkTrue(all(c("compounds","descriptor_types","descriptors") %in% dbListTables(conn)))
	checkTrue(inherits(conn,"DBIConnection"))
	dbDisconnect(conn)
}

test_ba.loadSdf<-function(){
	#DEACTIVATED("temp")
	data(sdfsample)


	conn = initDb("test1.db")
#	print("loading first half, no features, with exception")
#	checkException(loadSdf(conn,sdfsample[c(1,2,1,3)]))
	compIds=loadSdf(conn,sdfsample[c(1,2,3)])
	checkEquals(length(compIds),3)
	dbDisconnect(conn)


	conn = initDb("test1.db") # use new conn to make sure changes are durable
	print("loading first half, no features")
	compoundCount = dbGetQuery(conn,"SELECT count(*) FROM
										compounds WHERE format!='junk'")[1][[1]]
	checkEquals(compoundCount, 3)


	firstHalf=tempfile()
	write.SDF(sdfsample[1:50],firstHalf)
	secondHalf=tempfile()
	write.SDF(sdfsample[51:100],secondHalf)

	conn = initDb("test2.db")

	print("loading first half,with features")
	loadSdf(conn,firstHalf,function(sdfset)
			  data.frame( 
					  MW=MW(sdfset),
					  rings(sdfset, type="count", upper=6, arom=TRUE)) )

	print("loading incomplete features")
	checkException(loadSdf(conn,secondHalf,function(sdfset)
			  data.frame( CID=cid(sdfset),
					  rings(sdfset, type="count", upper=6, arom=TRUE)) ))

	print("loading second half")
	loadSdf(conn,secondHalf,function(sdfset)
			  data.frame( CID=cid(sdfset),
					  MW=MW(sdfset),
					  rings(sdfset, type="count", upper=6, arom=TRUE)) )
	print("done loading")


	compoundCount = dbGetQuery(conn,"SELECT count(*) FROM
										compounds WHERE format!='junk'")[1][[1]]
	checkEquals(compoundCount ,length(cid(sdfsample)))
	featureCount= dbGetQuery(conn,"SELECT count(*) FROM feature_mw")[1][[1]]
	checkEquals(featureCount ,length(cid(sdfsample)))
	dbDisconnect(conn)

	#test loading descriptors
	conn=initDb("test_desc.db")
	loadSdf(conn,firstHalf,function(sdfset) data.frame(MW=MW(sdfset)),
			  descriptors=function(sdfset) 
				data.frame(descriptor_type="ap",descriptor=unlist(lapply(ap(sdf2ap(sdfset)),
														function(ap) paste(ap,collapse=", ")))))
	descriptorCount= dbGetQuery(conn,"SELECT count(*) FROM descriptors")[1][[1]]
	checkEquals(descriptorCount,50)
	typeCount= dbGetQuery(conn,"SELECT count(*) FROM descriptor_types WHERE descriptor_type = 'ap' ")[1][[1]]
	checkEquals(typeCount,1)
}
test_bn.addNewFeatures<-function(){

	conn = initDb("test2.db")
	addNewFeatures(conn,function(sdfset) data.frame(new=cid(sdfset),new2=cid(sdfset)))
	features = listFeatures(conn)
	checkTrue("new" %in% features)
	checkTrue("new2" %in% features)

	#tables = dbListTables(conn)
	#checkTrue("feature_new" %in% tables)
	#checkTrue("feature_new2" %in% tables)
}

test_ca.findCompounds<-function(){
	#DEACTIVATED("temp")
	conn = initDb("test2.db")

	indexes = findCompounds(conn,"MW",c("MW < 400"))
	print(paste("found",length(indexes)," compounds"))
	checkEquals(length(indexes),70)

	checkException(findCompounds(conn,"MW",c("mw < 400","rings > 3")))

	indexes=findCompounds(conn,c("MW","RINGS"),c("mw < 400","rings > 3"))
	print(paste("found",length(indexes)," compounds"))
	checkEquals(length(indexes),20)

	dbDisconnect(conn)
}

test_da.getCompounds<-function(){
	#DEACTIVATED("temp")
	conn = initDb("test2.db")

	indexes = findCompounds(conn,"MW","mw < 400")

	sdfset = getCompounds(conn,indexes,keepOrder=TRUE)
	checkEquals(length(cid(sdfset)),70)
	ids2=findCompoundsByName(conn,sdfid(sdfset),keepOrder=TRUE)
	names(ids2)=c()
	checkEquals(indexes,ids2)
	
	getCompounds(conn,indexes,file="test.sdf")
	checkTrue(file.exists("test.sdf"))
	sdfFromFile = read.SDFset("test.sdf")
	checkEquals(length(cid(sdfFromFile)),70)
	dbDisconnect(conn)
}



test_ea.comparison <- function() {
	DEACTIVATED("local test")
	#filename = "~/runs/kinase/kinase.sdf"
	filename = "~/runs/protein/proteins.sdf"
	#filename = "~/runs/protein/proteins-1000.sdf"
#	options(warn=2)
	options(error=traceback)
	streamTest <- function(){
		#sink("/dev/null")
		print(system.time(sdfStream(input=filename,output="index.stream", silent=TRUE, fct=function(sdfset)
					 cbind(MW=MW(sdfset)  ))))
		print(system.time(index <- read.delim("index.stream",row.names=1)))
		#index[index$sdfid %in% c("3540","5329468","32014"),]
		print(system.time(queryIds <- index[ !is.na(index$MW) & index$MW < 400,]))
		#print(system.time(read.SDFindex(file=filename,index=queryIds,type="file",outfile="stream_result.sdf")))
		#print(system.time(sdfset<<-read.SDFindex(file=filename,index=queryIds)))
		#print(length(sdfset))

		#sink()
		#print(t1) #loadSdf
		#print(t2) #
		#print(t3) #findCompounds
		#print(t4) #write compounds
	}
	dbTest <- function(){
		print(system.time(conn<-initDb("tempdb")))
		print(system.time(loadSdf(conn,filename,function(sdfset)cbind(MW=MW(sdfset)),
										  function(x) {
											  #print(paste("apset time",length(x)))
											  #print(system.time(aps<<-sdf2ap(x)))
											  #data.frame(descriptor_type="ap",
															 #descriptor = unlist(lapply(ap(aps), 
																					#function(x) paste(x,collapse=", "))))
											  data.frame(descriptor_type=c(),descriptor=c())
										  }
										)))
		print(system.time(indexes <<- findCompounds(conn,"MW","mw < 400")))
		print(system.time(getCompounds(conn, indexes,file="dbtest_result.sdf")))
		print(system.time(sdfset<<-getCompounds(conn, indexes)))
		print(length(sdfset))
		dbDisconnect(conn)
	}

	unlink("tempdb")

	print("stream:")
	#Rprof()
	print(system.time(streamTest()))
	#Rprof(NULL)
	#summaryRprof("Rprof.out")
	
	print("db:")
#	Rprof()
	print(system.time(dbTest()))
}

test_fa.parBatchByIndex <- function(){
	message("starting parBatchByIndex")
	require(snow)
	cl = makeSOCKcluster(3)
	conn = initDb("test2.db")
	ids = dbGetQuery(conn,"SELECT compound_id FROM compounds WHERE format!='junk'")$compound_id
	print(str(ids))
	outfile="parBatch.out"

	unlink(outfile)
	unlink("parBench-*")
	parBatchByIndex(ids,batchSize = 10,cl=cl,
			   indexProcessor = function(indexSet,jobId){
					filename = paste("parBench-sub",jobId,sep="-")
					cat(indexSet,"\n", sep=" ",file=filename)
					filename
				},
				reduce =function(results) {

					print(paste("results: ",paste(results,collapse=",")))
					for(filename in results)
						cat(scan(filename,quiet=TRUE),sep="\n",file=outfile,append=TRUE)
				})
	resultIds = scan(outfile,quiet=TRUE)
	checkEquals(ids,resultIds)
}

test_ga.addDups <- function() {

	data(sdfsample)
	conn = initDb("test1.db")
	print("loading duplications")

	descFn = function(sdfset) data.frame(descriptor = paste("descriptor for ",sdfid(sdfset)),
													 descriptor_type = "testing")


	#add 3 dups by checksum
	count1= getCompoundCount(conn)
	loadSdf(conn,sdfsample[c(1,2,3)],descriptors=descFn)
	count2= getCompoundCount(conn)
	checkEquals(count1,count2)

	compIds = findCompoundsByName(conn,sdfid(sdfsample[1:3]),keepOrder=TRUE,allowMissing=TRUE)
	checkEquals(length(compIds),3)


	# add two dups and one update by checksum
	atomblock(sdfsample)[[1]][,]=8
	loadSdf(conn,sdfsample[c(1,2,3)],descriptors=descFn)
	count2= getCompoundCount(conn)
	checkEquals(count1+1,count2)
	compIds = findCompoundsByName(conn,sdfid(sdfsample[1]),allowMissing=TRUE)
	checkEquals(length(compIds),2)

	# add one dup and update by name
	atomblock(sdfsample)[[2]][,]=8
	loadSdf(conn,sdfsample[c(2,3)],descriptors=descFn, updateByName=TRUE)
	count2= getCompoundCount(conn)
	checkEquals(count1+1,count2)
}
getCompoundCount  <- function(conn){

	dbGetQuery(conn,"SELECT count(*) FROM compounds WHERE format!='junk'")[1][[1]]
}

test_ea.dupDescriptors <- function() {

	DEACTIVATED("local test only")
	#   descriptor_id          compound_id
	# 13881056,13881075      13884414,13884433
	# 41764092,41764082      41780498,41780494
   print("loading  duplicated descriptors")
	conn = initDb("dup_test.db")
	sdfs = read.SDFset("/home/khoran/ChemmineR/inst/unitTests/descriptor_dups.sdf")
	loadSdf(conn,sdfs,function(sdfset) data.frame(MW=MW(sdfset)),
			  descriptors=function(sdfset) 
				data.frame(descriptor_type="ap",descriptor=unlist(lapply(ap(sdf2ap(sdfset)),
														function(ap) paste(ap,collapse=", ")))))
	print("done loading")
	descriptorCount= dbGetQuery(conn,"SELECT count(*) FROM descriptors")[1][[1]]
	checkEquals(descriptorCount,1)
	typeCount= dbGetQuery(conn,"SELECT count(*) FROM descriptor_types WHERE descriptor_type = 'ap' ")[1][[1]]
	checkEquals(typeCount,1)
	linkCount= dbGetQuery(conn,"SELECT count(*) FROM compound_descriptors")[1][[1]]
	checkEquals(linkCount,7)
   dbDisconnect(conn)
	
}
