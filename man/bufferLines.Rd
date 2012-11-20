\name{bufferLines}
\alias{bufferLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Buffer File Input
}
\description{
	Buffer the input of files to increase efficiency
}
\usage{
bufferLines(fh, batchSize, lineProcessor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fh}{
%%     ~~Describe \code{fh} here~~
}
  \item{batchSize}{
%%     ~~Describe \code{batchSize} here~~
}
  \item{lineProcessor}{
%%     ~~Describe \code{lineProcessor} here~~
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
function (fh, batchSize, lineProcessor) 
{
    while (TRUE) {
        lines = readLines(fh, n = batchSize)
        if (length(lines) > 0) 
            lineProcessor(lines)
        else break
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
