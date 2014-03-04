\name{write.SMI}
\alias{write.SMI}
\title{
SMI export function
}
\description{
Writes one or many molecules stored in a \code{SMIset} object to a SMILES file.
}
\usage{
write.SMI(smi, file, cid = TRUE, ...)
}
\arguments{
  \item{smi}{
object of class \code{SMIset}
}
  \item{file}{
name of SMILES file to write to
}
  \item{cid}{
if \code{cid = TRUE} the compound identifiers will be exported by appending them in tab-separated format to each SMILES string
}
  \item{\dots}{
option to pass on additional arguments
}
}
\details{
...
}
\references{
SMILES (Simplified molecular-input line-entry system) format definition:
http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{write.SDF}
}
\examples{
## Instance of SMIset class
data(smisample); smiset <- smisample

## Write objects of classes SMIset to file with and 
## without compound identifiers
# write.SMI(smiset[1:4], file="sub.smi", cid=TRUE)
# write.SMI(smiset[1:4], file="sub.smi", cid=FALSE)

}
\keyword{ utilities }
