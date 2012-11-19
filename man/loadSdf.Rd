\name{loadSdf}
\alias{loadSdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Load SDF Data
}
\description{
	Load an SDF formatted file into the database
}
\usage{
loadSdf(conn, sdfFile, fct = function(x) data.frame(), descriptors=function(x) data.frame(descriptor=c(),descriptor_type=c()),
				Nlines = 10000, startline = 1, restartNlines = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{
	}
  \item{sdfFile}{
	}
  \item{fct}{
	}
	\item{descriptors}{
	}
  \item{Nlines}{
	}
  \item{startline}{
  }
  \item{restartNlines}{
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
function (conn, sdfFile, batchSize = 10000, validate = FALSE) 
{
    f = file(sdfFile, "r")
    tryCatch({
        dbOp(dbGetQuery(conn, "BEGIN TRANSACTION"))
        compoundLines = rep("", batchSize)
        compoundQueue = data.frame(name = rep(NA, 1000), definition = NA, 
            format = NA)
        compoundCount = 1
        lineNum = 1
        bufferLines(f, batchSize, function(lines) {
            for (line in lines) {
                compoundLines[lineNum] <<- line
                if (lineNum >= batchSize) 
                  compoundLines[lineNum * 2] <<- NA
                if (line == "$$$$") {
                  tryCatch({
                    if (validate) 
                      read.SDFset(compoundLines[1:lineNum])
                    def = paste(compoundLines[1:lineNum], collapse = "\n")
                    compoundQueue[compoundCount, ] <<- c(compoundLines[1], 
                      def, "sdf")
                    compoundCount <<- compoundCount + 1
                    if (compoundCount > dim(compoundQueue)[1]) {
                      print("loading batch")
                      print(compoundCount)
                      loadDb(conn, compoundQueue)
                      compoundCount <<- 1
                    }
                  }, error = function(e) {
                    print(paste("bad def found:", e$message))
                  })
                  lineNum <<- 0
                }
                lineNum <<- lineNum + 1
            }
        })
        print(paste("loading last batch", compoundCount))
        loadDb(conn, compoundQueue[1:compoundCount - 1, ])
        dbOp(dbCommit(conn))
    }, error = function(e) {
        print(paste("import failed:", e$message))
        traceback()
        dbOp(dbRollback(conn))
    })
    close(f)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
