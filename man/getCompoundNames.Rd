\name{getCompoundNames}
\alias{getCompoundNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Get Compound Names
}
\description{
	Fetch the names of the given compound ids, if they exist
}
\usage{
getCompoundNames(conn, compoundIds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
%%     ~~Describe \code{conn} here~~
}
  \item{compoundIds}{
%%     ~~Describe \code{compoundIds} here~~
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
function (conn, compoundIds) 
{
    result = selectInBatches(conn, compoundIds, function(ids) paste("SELECT compound_id, name FROM compounds where compound_id in (", 
        paste(ids, collapse = ","), ")"))
    as.matrix(merge(data.frame(compound_id = compoundIds), result)[[2]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
