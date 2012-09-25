\name{sdfStream}
\alias{sdfStream}
\title{
Streaming through large SD files
}
\description{
Streaming function to compute descriptors for large SD Files without consuming much memory. In addition to descriptor values, it returns a line index that defines the positions of each molecule in the source SD File. This line index can be used by the \code{read.SDFindex} function to retrieve specific compounds of interest from large SD Files without reading the entire file into memory.
}
\usage{
sdfStream(input, output, fct, Nlines = 10000, silent = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
file name of input SD file 
}
  \item{output}{
file name of tabular descriptor file
}
  \item{fct}{
function to select descriptor sets; any combination of descriptors, supported by \code{ChemmineR}, can be chosen here, as long as they can be represented in tabular format.
}
  \item{Nlines}{
number of lines to read from input SD File at a time; the memory consumption will be proportional to this value.
}
  \item{silent}{
if \code{silent=FALSE}, the processing status will be printed to the screen, while \code{silent=TRUE} suppresses this output.
}
  \item{\dots}{
Arguments to be passed to/from other methods.
}
}
\details{
...
}
\value{
Writes a descriptor matrix to a tabular file. The first and last line number (position index) of each molecule is specified in the first two columns of the tabular output file, respectively.
}
\references{
SDF format definition: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}
\seealso{
Import/export functions: \code{read.AP}, \code{read.SDFset}, \code{read.SDFstr}, \code{read.SDFstr}, \code{read.SDFset}, \code{write.SDFsplit}
}
\examples{
## Load sample data
library(ChemmineR)
data(sdfsample); sdfset <- sdfsample
\dontrun{write.SDF(sdfset, "test.sdf")

## Define descriptor set in a simple function
desc <- function(sdfset) {
        cbind(SDFID=sdfid(sdfset), 
              # datablock2ma(datablocklist=datablock(sdfset)), 
              MW=MW(sdfset), 
              groups(sdfset), 
              # AP=sdf2ap(sdfset, type="character"),
              rings(sdfset, type="count", upper=6, arom=TRUE)
        )
}

## Run sdfStream with desc function and write results to a file called 'matrix.xls'
sdfStream(input="test.sdf", output="matrix.xls", fct=desc, Nlines=1000)

## Select molecules from SD File using line index from sdfStream
indexDF <- read.delim("matrix.xls", row.names=1)[,1:4]
indexDFsub <- indexDF[indexDF$MW < 400, ] # Selects molecules with MW < 400
sdfset <- read.SDFindex(file="test.sdf", index=indexDFsub, type="SDFset")

## Write result directly to SD file without storing larger numbers of molecules in memory
read.SDFindex(file="test.sdf", index=indexDFsub, type="file", outfile="sub.sdf")

## Read atom pair string representation from file into APset
apset <- read.AP(file="matrix.xls", colid="AP")
cid(apsdf) <- as.character(indexDF$SDFID)  
}
}
\keyword{ utilities }
