\name{batchByIndex}
\alias{batchByIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Batch by Index
}
\description{
	When doing a select were the condition is a large number of ids it is not always possible
	to include them in a single SQL statement. This function will break the list of ids into chunks 
	and allow the indexProcessor to deal with just a small number of ids.
}
\usage{
batchByIndex(allIndices, indexProcessor, batchSize = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allIndices}{
		A vector of values that will be broken into batches and passed as an argument to the
		\code{indexProcessor} function.
	}
  \item{indexProcessor}{
  		A function that takes one batch if indices. It is called once for each batch.
		The return value from this function is ignored. To accumulate results you can	
		write to a global variable using the "<<-" operator. 
	}
  \item{batchSize}{
		The size of each batch. The last batch may be smaller than this value.
	}
}
\value{
	No value is returned.
}
\author{
	Kevin Horan
}

\seealso{
\code{\link{parBatchByIndex}}
}
\examples{
	
	\dontrun{
		result=NA
		indices = 1:10000
		
		#run a query on each batch of indexes, appending each result to
		# "result" as we go.
		batchByIndex(indices, function(indexBatch){
				df = dbGetQuery(dbConnection, generateQuery(indexBatch))
				result <<- if(is.na(result)) df else  rbind(result,df)
		},1000)
	}


}
