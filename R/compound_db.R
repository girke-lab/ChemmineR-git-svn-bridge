
#library(DBI)
#library(RSQLite)
#library(ChemmineR)

ddl <- function(conn,ddlSql){
#	tryCatch({
		dbGetQuery(conn,ddlSql)
		#rs=dbSendQuery(conn,ddlSql)
		#dbClearResult(rs)

#		},error=function(e){
#			if(length(grep("column definition is not unique",e$message))==0)
#				stop(paste("error sending query:",e$message))
#		})
}

initDb <- function(handle){

	if(is.character(handle)){ #assume sqlite filename
		require(RSQLite)
		driver = dbDriver("SQLite")
		conn = dbConnect(driver,dbname=handle)
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

		Map(function(sql) ddl(conn,sql),statements)
	}
	conn
}

loadDb <- function(conn,definitions,format,names=NA){

	if(inherits(conn,"SQLiteConnection")){
		insertDef <- function(defs)  dbGetPreparedQuery(conn,"INSERT INTO compounds(definition,format) VALUES(?,?)",
																		bind.data=data.frame(definition=defs,format="sdf"))
		insertNamedDef <- function(names,defs)  dbGetPreparedQuery(conn,"INSERT INTO compounds(name,definition,format) VALUES(?,?,?)",
																		bind.data=data.frame(name=names,definition=defs,format="sdf"))
	}else{
		insertDef <- function(defs)
			Map(function(def) dbGetQuery(conn, paste("INSERT INTO compounds(definition,format) VALUES('",def,"','",format,"')",
																  sep="")),defs)
		insertNamedDef <- function(names,defs)
			Map(function(def) dbGetQuery(conn, paste("INSERT INTO compounds(name,definition,format) VALUES('",
																	  name,"','",def,"','",format,"')",
																  sep="")),names,defs)
	}


	tryCatch({
		if(length(names)==1 &&  is.na(names) ){
			#if(length(definitions)==1)
				insertDef(definitions)
			#else
				#Map(insertDef, definitions)
		}else {
			if(length(names) != length(definitions))
				stop("names do not line up with definitions. 
					  Names can be blank, but there must be a string for each definition")

			#if(length(definitions)==1) #avoid some overhead of Map
				insertNamedDef(names,definitions)
			#else
				#Map(insertNamedDef, names,definitions)
		}
	},error=function(e){
			if(length(grep("column definition is not unique",e$message))==0)
				stop(paste("error sending query:",e$message))
	})
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
		chunk = fetch(rs,n=1000)
		if(length(chunk)==0)
			break;
		apply(chunk,1,rsProcessor)
	}
}
definition2SDFset <- function(def){
	#print(paste("def: ",def))
	#print(unlist(strsplit(def,"\n"),use.names=FALSE))
	#print("----------------------------------------------------------------")
	x = if(length(def) != 1) def else unlist(strsplit(def,"\n",fixed=TRUE),use.names=FALSE)
	read.SDFset(x)
}


loadSdf <- function(conn,sdfFile,batchSize=10000){
	f = file(sdfFile,"r")

	tryCatch({
		dbGetQuery(conn,"BEGIN TRANSACTION")
		compoundLines=rep("",batchSize)
		lineNum=1
		bufferLines(f,batchSize,function(lines){
					for(line in lines){
						#print(paste(lineNum,":",line))
						compoundLines[lineNum]<<-line
						if(lineNum >= batchSize)
							compoundLines[lineNum*2]<<-NA #expand array

						if(line == "$$$$"){ #end of a compound
							name=compoundLines[1]

							def=paste(compoundLines[1:lineNum],collapse="\n")
							#tryCatch(definition2SDFset(compoundLines[1:lineNum]),
							tryCatch(definition2SDFset(def),
										error=function(e) {
											print(paste("bad def found:",e$message))
											print(def)
											print("context:")
											print(compoundLines)
										})

							loadDb(conn,def,"sdf",name)
							lineNum<<-0
						}
						lineNum <<- lineNum + 1
					}
				 })
		dbCommit(conn)
	},error=function(e){
		print(paste("import failed:",e$message))
		dbRollback(conn)
	})
	close(f)
}

loadSmiles <- function(conn, smileFile,batchSize=10000){

	f = file(smileFile,"r")
	bufferLines(f,batchSize=batchSize,function(lines) loadDb(conn,lines,"smile"))
	close(f)
}
findCompounds <- function(conn,test){
	rs = dbSendQuery(conn,"SELECT compound_id,definition FROM compounds ")
	matches = c()
	bufferResultSet(rs,function(row){
				tryCatch({
						if(test(definition2SDFset(row[2])[[1]])) #if definition passes test
							matches <<- c(matches,row[1]) #record id number
					},error=function(e) print(paste("error:",e$message,"on",row[1]))
				)
			 })
	dbClearResult(rs)
	as.numeric(matches)
}
getCompounds <- function(conn,compoundIds,filename=NA){
	
	if(!is.na(filename))
		f=file(filename,"w")

	sdfset = c()

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
		bufferResultSet(rs,function(row){
					sdf = definition2SDFset(row[2])
					cid(sdf) = row[1]
					if(!is.na(filename))
						write.SDF(sdf,file=f)
					else
						sdfset <<- if(length(sdfset)==0)
									sdf
								else
									c(sdfset,sdf)
				 })
		dbClearResult(rs)
	}


	if(!is.na(filename))
		close(f)
	sdfset
}


