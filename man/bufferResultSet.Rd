\name{bufferResultSet}
\alias{bufferResultSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Buffer Query Results
}
\description{
	Allow query results to be processed in batches for efficiency.
}
\usage{
bufferResultSet(rs, rsProcessor, batchSize = 1000,closeRS=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{
		A DBIResult object, usually from \code{dbSendQuery}.
	}
  \item{rsProcessor}{
	  Each batch will be passed as a data frame to this function for processing.
	}
  \item{batchSize}{
		The number of rows to read in each batch
	}
	\item{closeRS}{
		Should the result set be closed by this function when it is done?
	}
}
\value{
	No value.
}
\author{
	Kevin Horan
}
\examples{
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rs, rsProcessor, batchSize = 1000) 
{
    while (TRUE) {
        chunk = fetch(rs, n = batchSize)
        if (dim(chunk)[1] == 0) 
            break
        rsProcessor(chunk)
    }
  }
}
