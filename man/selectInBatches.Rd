\name{selectInBatches}
\alias{selectInBatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Select in Batches
}
\description{
	When doing a select were the condition is a large number of ids it is not always possible
	to include them in a single SQL statement. This function will break the list of ids into chunks 
	and  send the query for each batch. The resutls are appended and returned as one data frame.
}
\usage{
selectInBatches(conn, allIndices, genQuery, batchSize = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
		Database connection object
	}
  \item{allIndices}{
	  A vector of indices to pass to the genQuery function in batches.
	}
  \item{genQuery}{
	  A function which takes a vector of indices and constructs an SQL SELECT statement
	  returning records for the given indicies.
	}
  \item{batchSize}{
	  How many indicies to put in each batch.
	}
}
\value{
	A data frame with the results of the query as if all inidices had been included
	in a single SELEcT statement.
}
\author{
	Kevin Horan
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (conn, allIndices, genQuery, batchSize = 1e+05) 
{
    batchByIndex(allIndices, function(indexBatch) {
        df = dbGetQuery(conn, genQuery(indexBatch))
        result = rbind(result, df)
    }, batchSize)
    result
  }
}
