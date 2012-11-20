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
%%     ~~Describe \code{allIndices} here~~
}
  \item{indexProcessor}{
%%     ~~Describe \code{indexProcessor} here~~
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
function (allIndices, indexProcessor, batchSize = 1e+05) 
{
    numIndices = length(allIndices)
    if (numIndices == 0) 
        return()
    start = 1
    for (end in seq(1, numIndices, by = batchSize) + batchSize) {
        end = min(end, numIndices)
        if (debug) 
            print(paste(start, end))
        indexSet = allIndices[start:end]
        start = end + 1
        indexProcessor(indexSet)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
