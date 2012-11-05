\name{loadSmiles}
\alias{loadSmiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Load SMILES Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadSmiles(conn, smileFile, batchSize = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
%%     ~~Describe \code{conn} here~~
}
  \item{smileFile}{
%%     ~~Describe \code{smileFile} here~~
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
function (conn, smileFile, batchSize = 10000) 
{
    f = file(smileFile, "r")
    bufferLines(f, batchSize = batchSize, function(lines) loadDb(conn, 
        data.frame(definition = lines, format = "smile")))
    close(f)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
