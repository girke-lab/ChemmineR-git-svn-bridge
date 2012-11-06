
#library(DBI)
#library(RSQLite)
#library(ChemmineR)


dbOp<-function(dbExpr){
	#print(as.character(substitute(dbExpr)))
	#print(system.time(dbExpr))
	dbExpr
}

initDb <- function(handle){

	if(is.character(handle)){ #assume sqlite filename
		require(RSQLite)
		driver = dbDriver("SQLite")
		conn = dbConnect(driver,dbname=handle,cache_size=100000)
	}else if(inherits(handle,"DBIConnection")){ #the user passed in a connection object
		conn=handle
	}else{
		stop("handle must be a SQLite database name or a DBIConnection")
	}

	tableList=dbListTables(conn)
	#print(tableList)

	if( ! all(c("compounds","descriptor_types","descriptors") %in% tableList)) {
		print("createing db")

		statements = unlist(strsplit(paste(readLines("inst/schema/compounds.SQLite"),collapse=""),";",fixed=TRUE))
		#print(statements)

		Map(function(sql) dbOp(dbGetQuery(conn,sql)),statements)
	}
	conn
}

#loadDb <- function(conn,definitions,format,names=NA){
loadDb <- function(conn,data,featureGenerator){

	print(paste("loading ",paste(dim(data),collapse=" "),"compounds"))

	data=cbind(data,definition_checksum=sapply(as.vector(data[,"definition"]),
													  function(def) digest(def,serialize=FALSE) ))

	tryCatch({
		addNewFeatures(conn,data,featureGenerator)

		if(all(c("name","definition","format") %in% colnames(data))){
			insertNamedDef(conn,data)
		}else if(all(c("definition","format") %in% colnames(data))){
			insertDef(conn,data)
		}else {
			stop("given data must have columns either (definition,format), or (name,definiton,format)")
		}

		##updateFeatures(conn,data,featureGenerator)
		insertUserFeatures(conn,data)
	},error=function(e){
			if(length(grep("column definition is not unique",e$message))==0){
				#print(sys.calls())
				stop(paste("error sending query:",e$message))
			}
	})
}
addNewFeatures <- function(conn,data, featureGenerator){
	if(dim(data)[1] == 0){ # no data given
		warning("no data given to addNewFeatures, doing nothing")
		return()
	}

	compoundFields = dbListFields(conn,"compounds")
	userFieldNames = setdiff(colnames(data),compoundFields)

	tableList=dbListTables(conn)
	existingFeatures = sub("^feature_","",tableList[grep("^feature_.*",tableList)])

	missingFeatures = setdiff(existingFeatures,userFieldNames)
	if(length(missingFeatures) != 0)
		stop(paste("missing features:",paste(missingFeatures,collapse=",")))

	newFeatures = setdiff(userFieldNames,existingFeatures)

	#will need to query all existing compounds and add these features
	if(length(newFeatures) != 0){
		message("adding new features to existing compounds. This could take a while")
		lapply(newFeatures,function(name) createFeature(conn,name,class(data[,name])))
		indexExistingCompounds(conn,newFeatures,featureGenerator)
		
	}


}
insertUserFeatures<- function(conn,data){
	print("inserting user features")

	if(dim(data)[1] == 0){ # no data given
		warning("no data given to insertUserFeatures, doing nothing")
		return()
	}


	compoundFields = dbListFields(conn,"compounds")
	userFieldNames = setdiff(colnames(data),compoundFields)

	tableList=dbListTables(conn)
	existingFeatures = sub("^feature_","",tableList[grep("^feature_.*",tableList)])

		
	#index all new compounds for all features

	#fetch newly inserted compounds
	
	df = dbGetQuery(conn,paste("SELECT compound_id,definition_checksum FROM compounds LEFT JOIN feature_",
								 userFieldNames[1]," as f USING(compound_id)  WHERE f.compound_id IS
								 NULL",sep=""))
	data= merge(data,df)
	sapply(userFieldNames,function(name) insertFeature(conn,name,data))





}
bufferLines <- function(fh,batchSize,lineProcessor){
	while(TRUE){
		lines = readLines(fh,n=batchSize)
		if(length(lines) > 0)
			lineProcessor(lines)
		else
			break;
	}
}
bufferResultSet <- function(rs,rsProcessor,batchSize=1000){
	while(TRUE){
		chunk = fetch(rs,n=batchSize)
		if(dim(chunk)[1]==0) # 0 rows in data frame
			break;
		#apply(chunk,1,rsProcessor)
		rsProcessor(chunk)
	}
}
definition2SDFset <- function(defs){

#	f = file()
#	lapply(defs,function(def) cat(def,file=f))
#	flush(f)
#	sdfset = read.SDFset(read.SDFstr(f))
#	close(f)
#	return(sdfset)

	read.SDFset(unlist(strsplit(defs,"\n",fixed=TRUE)))


#	sdfset=c()
#	count=0
#	suppressWarnings(lapply(defs,function(def){ 
#			 sdf=read.SDFset(unlist(strsplit(def,"\n",fixed=TRUE),use.names=FALSE))
#			 #cid(sdf)=sdfid(sdf)
#			 #print(paste(cid(sdf),sdfid(sdf)))
#			 sdfset<<-if(length(sdfset)==0) sdf else c(sdfset,sdf)
#			 count<<-count+1 #this somehow keeps the memory under control. who knows why
#			# l=length(sdfset)
#			# if( l %% 500 == 0){
#			#	 print(length(sdfset))
#			#	 gc()
#			# }
#	} ))
#	sdfset
}

loadSdf2 <- function(conn,sdfFile,fct=function(x) cbind(), Nlines=10000, startline=1, restartNlines=100000){
	## Define loop parameters 
	stop <- FALSE 
	f <- file(sdfFile, "r")
	n <- Nlines
	offset <- 0
	## For restarting sdfStream at specific line assigned to startline argument. If assigned
        ## startline value does not match the first line of a molecule in the SD file then it 
        ## will be reset to the start position of the next molecule in the SD file.
	if(startline!=1) { 
		fmap <- file(sdfFile, "r")
		shiftback <- 2
		chunkmap <- scan(fmap, skip=startline-shiftback, nlines=restartNlines, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n")
		startline <- startline + (which(grepl("^\\${4,4}", chunkmap, perl=TRUE))[1] + 1 - shiftback)
		if(is.na(startline)) stop("Invalid value assigned to startline.")
		dummy <- scan(f, skip=startline-2, nlines=1, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n")
		close(fmap)
		offset <- startline - 1 # Maintains abolut line positions in index
	}
	counter <- 0
	cmpid <- 1
	partial <- NULL
	dbOp(dbGetQuery(conn,"BEGIN TRANSACTION"))
	while(!stop) {
		counter <- counter + 1
		chunk <- scan(f, n=n, what="a", blank.lines.skip=FALSE, quiet=TRUE, sep ="\n") # scan has more flexibilities for reading specific line ranges in files.
		if(length(chunk) > 0) {
			if(length(partial) > 0) {
				chunk <- c(partial, chunk)
			}
			## Assure that lines of least 2 complete molecules are stored in chunk if available
			inner <- sum(grepl("^\\${4,4}", chunk, perl=TRUE)) < 2
			while(inner) {
				chunklength <- length(chunk)
				chunk <- c(chunk, readLines(f, n = n))
				if(chunklength == length(chunk)) { 
					inner <- FALSE 
				} else {
					#inner <- sum(grepl("^\\${4,4}", chunk, perl=TRUE)) < 2
					inner <- sum(grepl("$$$$", chunk, fixed=TRUE)) < 2
				}
			}
			y <- regexpr("^\\${4,4}", chunk, perl=TRUE) # identifies all fields that start with a '$$$$' sign
			index <- which(y!=-1)
			indexDF <- data.frame(start=c(1, index[-length(index)]+1), end=index)
			complete <- chunk[1:index[length(index)]]
			if((index[length(index)]+1) <= length(chunk)) {
				partial <- chunk[(index[length(index)]+1):length(chunk)]
			} else {
				partial <- NULL
			}

			sdfset <- read.SDFset(read.SDFstr(complete))
			valid = validSDF(sdfset)
			sdfset=sdfset[valid]


			userFeatures <- fct(sdfset)

	#		sdfstrList=as(as(sdfset,"SDFstr"),"list")
	#		names = unlist(Map(function(x) x[1],sdfstrList))
	#		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )

			defs=apply(indexDF,1,function(row) paste(complete[row[1]:row[2]],collapse="\n"))
			names = complete[indexDF[,1]]
			systemFields=data.frame(name=names[valid],definition=defs[valid],format="sdf")

			allFields = if(length(userFeatures)!=0) cbind(systemFields,userFeatures) else systemFields
			loadDb(conn,allFields,fct)

			
		}
		if(length(chunk) == 0) {
			stop <- TRUE
			close(f)
		}
	}

	dbOp(dbCommit(conn))

}
loadSdf3 <- function(conn,sdfFile,batchSize=10000,validate=FALSE){

	dbOp(dbGetQuery(conn,"BEGIN TRANSACTION"))
	sdfStream(input=sdfFile,output="/dev/null",silent=TRUE,fct=function(sdfset){
		sdfstrList=as(as(sdfset,"SDFstr"),"list")
		names = unlist(Map(function(x) x[1],sdfstrList))
		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )
		data=data.frame(name=names,definition=defs,format="sdf")
		loadDb(conn,data)
		cbind(MW=1:length(sdfset))
	})
	dbOp(dbCommit(conn))
}

loadSdf <- function(conn,sdfFile,batchSize=10000,validate=FALSE){
	f = file(sdfFile,"r")

	tryCatch({
		dbOp(dbGetQuery(conn,"BEGIN TRANSACTION"))
		compoundLines=rep("",batchSize)
		compoundQueue=data.frame(name=rep(NA,1000),definition=NA,format=NA)
		compoundCount=1
		lineNum=1

		bufferLines(f,batchSize,function(lines){
					for(line in lines){
						#print(paste(lineNum,":",line))
						compoundLines[lineNum]<<-line
						if(lineNum >= batchSize)
							compoundLines[lineNum*2]<<-NA #expand array

						if(line == "$$$$"){ #end of a compound
							tryCatch({
								if(validate)
									read.SDFset(compoundLines[1:lineNum])

								def=paste(compoundLines[1:lineNum],collapse="\n")

								compoundQueue[compoundCount,]<<-c(compoundLines[1],def,"sdf")
								compoundCount<<- compoundCount+1
								if(compoundCount > dim(compoundQueue)[1]){
									print("loading batch")
									print(compoundCount)
									loadDb(conn,compoundQueue)
									compoundCount<<-1
								}
							},
							error=function(e) {
								print(paste("bad def found:",e$message))
								#print("context:")
								#print(compoundLines)
							})
							
							lineNum<<-0
						}
						lineNum <<- lineNum + 1
					}
				 })
		print(paste("loading last batch",compoundCount))
		#print(compoundQueue[,c(1,3)])
		loadDb(conn,compoundQueue[1:compoundCount-1,])
		dbOp(dbCommit(conn))
	},error=function(e){
		print(paste("import failed:",e$message))
		traceback()
		dbOp(dbRollback(conn))
	})
	close(f)
}

loadSmiles <- function(conn, smileFile,batchSize=10000){

	f = file(smileFile,"r")
	bufferLines(f,batchSize=batchSize,function(lines) loadDb(conn,data.frame(definition=lines,format="smile")))
	close(f)
}

findCompounds <- function(conn,featureNames,tests){

	# SELECT compound_id FROM compounds join f1 using(compound_id) join f2 using(compound_id)
	# ... where test1 AND test2 AND ...
	featureTables = paste("feature_",featureNames,sep="")
	sql = paste("SELECT compound_id FROM compounds JOIN ",paste(featureTables,collapse=" USING(compound_id) JOIN "),
					" USING(compound_id) WHERE ",paste("(",paste(tests,collapse=") AND ("),")") ) 
	
	#print(paste("query sql:",sql))
	result = dbGetQuery(conn,sql)
	result[1][[1]]

}
#findCompounds_slow <- function(conn,test){
#	rs = dbOp(dbSendQuery(conn,"SELECT compound_id,definition FROM compounds "))
#	matches = c()
#	bufferResultSet(rs,function(row){
#				tryCatch({
#						if(test(definition2SDFset(row[2])[[1]])) #if definition passes test
#							matches <<- c(matches,row[1]) #record id number
#					},error=function(e) print(paste("error:",e$message,"on",row[1]))
#				)
#			 },1)
#	dbOp(dbClearResult(rs))
#	as.numeric(matches)
#}
getCompounds <- function(conn,compoundIds,filename=NA){
	
	if(!is.na(filename)){
		f=file(filename,"w")
#	}else {
#		f=file()
#	}

		resultProcessor = function(rows){
			lapply(rows[2][[1]],function(def) cat(def,file=f))
			#cat(rows[2][[1]],file=f)
		}
	}else{
		sdfset = rep(NA,length(compoundIds))
		count=1
		resultProcessor = function(rows){
			print(dim(rows))
				#print(system.time(sdfs <<- definition2SDFset(rows[2][[1]])))
				#print(system.time(cid(sdfs) <<- as.character(rows[1][[1]])))
				#print(system.time(sdfset <<- if(length(sdfset)==0) sdfs else c(sdfset,sdfs)))
				print(system.time({
										l=length(rows[2][[1]])
										sdfset[count:(count+l-1)]<<-as(definition2SDFset(rows[2][[1]]),"SDF")
										names(sdfset)[count:(count+l-1)]<<-as.character(rows[1][[1]])
										count<<-count+l

				}))
		}
	}


	indexChunkSize=100000
	start=1
	print(paste("length:",length(compoundIds)))
	for(end in seq(1,length(compoundIds),by=indexChunkSize)+indexChunkSize){
	
		end = min(end,length(compoundIds))
		print(paste(start,end))
		compoundIdSet = compoundIds[start:end]
		start=end+1

		rs = dbSendQuery(conn,
							  paste("SELECT compound_id,definition FROM compounds where compound_id in (",
									  paste(compoundIdSet,collapse=","),")"))
		bufferResultSet(rs, resultProcessor,1000)
		dbClearResult(rs)
	}


	if(!is.na(filename)){
		close(f)
	}
	else{
	#	sdfset = read.SDFset(read.SDFstr(f))
		return(as(sdfset,"SDFset"))
	}
#	sdfset
}

indexExistingCompounds <- function(conn,newFeatures,featureGenerator){
	rs = dbOp(dbSendQuery(conn,"SELECT compound_id,definition FROM compounds "))
	bufferResultSet(rs,function(rows){
				tryCatch({
						sdfset = definition2SDFset(rows)
						userFeatures  = featureGenerator(sdfset)
						print(userFeatures)
						lapply(newFeatures, function(name){
								 print(paste("name:",name))
								 insertFeature(conn,name,cbind(compound_id=rows[1,],userFeatures[name]))
					   })
					},error=function(e) stop(paste("error in indexExistingCompounds:",e$message))
				)
			 })
	dbClearResult(rs)
}

createFeature <- function(conn,name, type){


	sqlType = if(type == "numeric") "NUMERIC" else "TEXT"
	print(paste("sql type: ",sqlType))
	dbGetQuery(conn,
		paste("CREATE TABLE feature_",name," (
			compound_id INTEGER REFERENCES compound(compound_id) ON DELETE CASCADE, ",
			name," ",sqlType," )",sep=""))
	print("made table")
	dbGetQuery(conn,paste("CREATE INDEX feature_",name,"_index ON
								 feature_",name,"(",name,")",sep=""))
	print("made index")

}

insertDef <- function(conn,data)  
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(definition,definition_checksum,format) ",
								 "VALUES(:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else{
		apply(data,1,function(row) dbOp(dbGetQuery(conn, 
						 paste("INSERT INTO compounds(definition,definition_checksum,format)
																			  VALUES('",row[1],"','",row[2],"','",row[3],"')", sep=""))))
	}

insertNamedDef <- function(conn,data) 
	if(inherits(conn,"SQLiteConnection")){
		print("preparing")
		dbGetPreparedQuery(conn,paste("INSERT INTO compounds(name,definition,definition_checksum,format) ",
								 "VALUES(:name,:definition,:definition_checksum,:format)",sep=""), bind.data=data)
	}else{
		apply(data,1,function(row) dbOp(dbGetQuery(conn, 
						 paste("INSERT INTO compounds(name,definition,definition_checksum,format) VALUES('",
																	  row[1],"','",row[2],"','",row[3],"','",row[4],"')", sep=""))))
	}

insertFeature <- function(conn,name,values){
	print(paste("name:",name))
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO feature_",name,"(compound_id,",name,") ",
												 "VALUES(:compound_id,:",name,")",sep=""), bind.data=values)
	}else{
		apply(data,1,function(row) 
				dbGetQuery(conn,paste("INSERT INTO feature_",name,"(compound_id,",name,")
											 VALUES(",row[1],",", if(!is.numeric(row[2]))
													  {paste("'",row[2],"'",sep="")} else {row[2]},")")))
	}
}
