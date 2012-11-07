\name{getCompounds}
\alias{getCompounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Get Compounds From Database
}
\description{
	Create SDF objects from the given set of compound ids. Id numbers can be found using the findCompounds function.
}
\usage{
getCompounds(conn, compoundIds, filename = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
%%     ~~Describe \code{conn} here~~
}
  \item{compoundIds}{
%%     ~~Describe \code{compoundIds} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
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
function (conn, compoundIds, filename = NA) 
{
    if (!is.na(filename)) 
        f = file(filename, "w")
    sdfset = c()
    indexChunkSize = 1e+05
    start = 1
    print(paste("length:", length(compoundIds)))
    for (end in seq(1, length(compoundIds), by = indexChunkSize) + 
        indexChunkSize) {
        end = min(end, length(compoundIds))
        print(paste(start, end))
        compoundIdSet = compoundIds[start:end]
        start = end + 1
        rs = dbOp(dbSendQuery(conn, paste("SELECT compound_id,definition FROM compounds where compound_id in (", 
            paste(compoundIdSet, collapse = ","), ")")))
        bufferResultSet(rs, function(row) {
            sdf = definition2SDFset(row[2])
            cid(sdf) = row[1]
            if (!is.na(filename)) 
                write.SDF(sdf, file = f)
            else sdfset <<- if (length(sdfset) == 0) 
                sdf
            else c(sdfset, sdf)
        }, 1)
        dbOp(dbClearResult(rs))
    }
    if (!is.na(filename)) 
        close(f)
    sdfset
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
