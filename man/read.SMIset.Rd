\name{read.SMIset}
\alias{read.SMIset}
\title{
SMILES file to \code{SMIset}
}
\description{
Imports one or many molecules from a SMILES file and stores content in a \code{SMIset} container. The input file is expected to
contain one SMILES string per row with tab-separated compound identifiers at the end of each line. The compound identifiers
are optional. 
}
\usage{
read.SMIset(file, removespaces = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
path/name to a SMILES file
}
  \item{removespaces}{
if set to \code{TRUE} spaces will be removed
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
...
}
\value{
\item{SMIset}{ for details see ?"SMIset-class"}
}
\references{
SMILES (Simplified molecular-input line-entry system) format definition:
http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{read.SDFset}
}
\examples{
## Write instance of SMIset class to SMILES file
data(smisample); smiset <- smisample
# write.SMI(smiset[1:4], file="sub.smi")

## Import SMILES file 
# read.SMIset("sub.smi")

}
\keyword{ utilities }
