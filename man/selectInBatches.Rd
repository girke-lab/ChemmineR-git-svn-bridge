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
%%     ~~Describe \code{conn} here~~
}
  \item{allIndices}{
%%     ~~Describe \code{allIndices} here~~
}
  \item{genQuery}{
%%     ~~Describe \code{genQuery} here~~
}
  \item{batchSize}{
%%     ~~Describe \code{batchSize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
