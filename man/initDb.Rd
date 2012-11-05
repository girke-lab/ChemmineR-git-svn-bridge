\name{initDb}
\alias{initDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Iinitialize Database
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initDb(handle)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{handle}{
%%     ~~Describe \code{handle} here~~
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
function (handle) 
{
    if (is.character(handle)) {
        require(RSQLite)
        driver = dbDriver("SQLite")
        conn = dbConnect(driver, dbname = handle, cache_size = 1e+05)
    }
    else if (inherits(handle, "DBIConnection")) {
        conn = handle
    }
    else {
        stop("handle must be a SQLite database name or a DBIConnection")
    }
    tableList = dbListTables(conn)
    if (!all(c("compounds", "descriptor_types", "descriptors") \%in\% 
        tableList)) {
        print("createing db")
        statements = unlist(strsplit(paste(readLines("inst/schema/compounds.SQLite"), 
            collapse = ""), ";", fixed = TRUE))
        Map(function(sql) dbOp(dbGetQuery(conn, sql)), statements)
    }
    conn
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
