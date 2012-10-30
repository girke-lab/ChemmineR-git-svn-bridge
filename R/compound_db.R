
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
loadDb <- function(conn,data){

	if(inherits(conn,"SQLiteConnection")){
		insertDef <- function(defs)  dbOp(dbGetPreparedQuery(conn,"INSERT INTO compounds(definition,format) VALUES(?,?)",
																		bind.data=data))
		insertNamedDef <- function(names,defs)  dbOp(dbGetPreparedQuery(conn,"INSERT INTO compounds(name,definition,format) VALUES(?,?,?)",
																		bind.data=data))
	}else{
		insertDef <- function(defs)
			apply(data,1,function(row) dbOp(dbGetQuery(conn, paste("INSERT INTO compounds(definition,format)
																			  VALUES('",row[1],"','",row[2],"')", sep=""))))
		insertNamedDef <- function(names,defs)
			apply(data,1,function(row) dbOp(dbGetQuery(conn, paste("INSERT INTO compounds(name,definition,format) VALUES('",
																	  row[1],"','",row[2],"','",row[3],"')", sep=""))))
	}

	#print("loading")
	#print(data)
	print(paste("loading ",paste(dim(data),collapse=" "),"compounds"))

	tryCatch({
		if(dim(data)[2]==2){
				insertDef(data)
		}else if(dim(data)[2]==3){
				insertNamedDef(data)
		}else {
			stop("given data must have columns either (definition,format), or (name,definiton,format)")
			#if(length(names) != length(definitions))
				#stop("names do not line up with definitions. 
					  #Names can be blank, but there must be a string for each definition")

			#if(length(definitions)==1) #avoid some overhead of Map
				#insertNamedDef(names,definitions)
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

loadSdf2 <- function(conn,sdfFile, Nlines=10000, startline=1, restartNlines=100000){
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
			valid <- validSDF(sdfset)
			#sdfset=sdfset[valid]

	#		sdfstrList=as(as(sdfset,"SDFstr"),"list")
	#		names = unlist(Map(function(x) x[1],sdfstrList))
	#		defs = unlist(Map(function(x) paste(x,collapse="\n"), sdfstrList) )

			defs=apply(indexDF,1,function(row) paste(complete[row[1]:row[2]],collapse="\n"))
			names = complete[indexDF[,1]]
			loadDb(conn,data.frame(name=names[valid],definition=defs[valid],format="sdf"))

			
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
									definition2SDFset(compoundLines[1:lineNum])

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
findCompounds <- function(conn,test){
	rs = dbOp(dbSendQuery(conn,"SELECT compound_id,definition FROM compounds "))
	matches = c()
	bufferResultSet(rs,function(row){
				tryCatch({
						if(test(definition2SDFset(row[2])[[1]])) #if definition passes test
							matches <<- c(matches,row[1]) #record id number
					},error=function(e) print(paste("error:",e$message,"on",row[1]))
				)
			 })
	dbOp(dbClearResult(rs))
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

		rs = dbOp(dbSendQuery(conn,
							  paste("SELECT compound_id,definition FROM compounds where compound_id in (",
									  paste(compoundIdSet,collapse=","),")")))
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
		dbOp(dbClearResult(rs))
	}


	if(!is.na(filename))
		close(f)
	sdfset
}


