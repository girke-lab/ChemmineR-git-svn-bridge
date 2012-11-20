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
%%     ~~Describe \code{rs} here~~
}
  \item{rsProcessor}{
%%     ~~Describe \code{rsProcessor} here~~
}
  \item{batchSize}{
%%     ~~Describe \code{batchSize} here~~
}
	\item{closeRS}{
		Should the result set be closed by this function when it is done?
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
