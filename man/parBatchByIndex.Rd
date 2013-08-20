\name{parBatchByIndex}
\alias{parBatchByIndex}
\title{
	Parallel Batch By Index
}
\description{
	Takes an index set, breaks it into batches and runs the given function on each batch
	in parallel using the given cluster. See \code{\link{batchByIndex}} for the non-parallel version.

	When doing a select were the condition is a large number of ids it is not always possible
	to include them in a single SQL statement. This function will break the list of ids into chunks 
	and allow the indexProcessor to deal with just a small number of ids.
}
\usage{
parBatchByIndex(allIndices, indexProcessor, reduce, cl, batchSize = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allIndices}{
		A vector of values that will be broken into batches and passed as an argument to the
		\code{indexProcessor} function.
	}
  \item{indexProcessor}{
  		A function that takes one batch if indices. It is called once for each batch, possibly in
		parallel. The return value of this function is collected into a list and passed to the
		\code{reduce} function after all jobs have finished.
	}
  \item{reduce}{
  		This function is run after all jobs have finished. It is called with a list of return values from
		the \code{indexProcessor} function runs.  The order of batchs is maintained. The return value of the \code{reduce}
		function is then returned.

		The idea is that this function merges all the results together into one result.
	}
  \item{cl}{
  		A SNOW cluster to run jobs on.
	}
  \item{batchSize}{
		The size of each batch. The last batch may be smaller than this value.
	}
}
\value{
	The return value of the \code{reduce} function is returned.
}
\author{
	Kevin Horan
}
\seealso{
\code{\link{batchByIndex}}
}
\examples{
	\dontrun{

		cl = makeCluster(2) # create a SNOW cluster

		#function to run a query for each batch of indexes
		job = function(indexBatch)
				dbGetQuery(dbConnection, paste("SELECT weight FROM table WHERE id IN (",paste(indexBatch,collapse=","),")"))

		# function to combine all the results, in this case by summing them up
		reduce = function(results) sum(unlist(results))

		indices = 1:10000

		#run queries in parallel and then sum the results
		totalWeight = parBatchByIndex(indices,job,reduce,cl, 1000)

	}
}
