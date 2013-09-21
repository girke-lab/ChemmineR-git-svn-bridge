\name{smisample}
\alias{smisample}
\docType{data}
\title{
SMILES file in \code{SMIset} object
}
\description{
First 100 compounds from PubChem SD file (Compound_00650001_00675000.sdf.gz) converted to SMILES format
}
\usage{data(smisample)}
\format{
  Object of class \code{smiset}
}
\details{
Object stores 100 molecules from a sample SMILES file.
}
\source{
ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_00650001_00675000.sdf.gz
}
\references{
SMILES (Simplified molecular-input line-entry system) format definition:
http://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
}
\examples{
data(smisample)
smiset <- smisample
view(smiset[1:4])
}
\keyword{datasets}
