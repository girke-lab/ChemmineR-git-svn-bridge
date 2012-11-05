\name{findCompounds}
\alias{findCompounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Find Compounds in Database
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findCompounds(conn, featureNames, tests)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
%%     ~~Describe \code{conn} here~~
}
  \item{featureNames}{
%%     ~~Describe \code{featureNames} here~~
}
  \item{tests}{
%%     ~~Describe \code{tests} here~~
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
function (conn, featureNames, tests) 
{
    featureTables = paste("feature_", featureNames, sep = "")
    sql = paste("SELECT compound_id FROM compounds JOIN ", paste(featureTables, 
        collapse = " USING(compound_id) JOIN "), " USING(compound_id) WHERE ", 
        paste("(", paste(tests, collapse = ") AND ("), ")"))
    result = dbGetQuery(conn, sql)
    result[1][[1]]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
