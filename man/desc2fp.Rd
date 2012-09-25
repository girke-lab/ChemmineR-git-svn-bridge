\name{desc2fp}
\alias{desc2fp}
\title{
Fingerprints from descriptor vectors
}
\description{
Generates fingerprints from descriptor vectors such as atom pairs. The obtained fingerprints can be used for structure similarity searching.
}
\usage{
desc2fp(x, descnames, type = "matrix")
}
\arguments{
  \item{x}{
Objects of classes \code{APset}
}
  \item{descnames}{
Names of descriptors (\code{character} vector) to consider for fingerprint encoding 
}
  \item{type}{
return fingerprint set as \code{matrix} or \code{character} vectors 
}
}
\details{
...
}
\value{
\code{matrix} or \code{character} vectors
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", J Chem Inf Comput Sci.
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{sdf2ap}, \code{SDF2apcmp}, \code{apset2descdb}, \code{cmp.search}, \code{cmp.similarity}

Related classes: SDF, SDFset, SDFstr, APset. 
}
\examples{
## Instance of SDFset class
data(sdfsample)
sdfset <- sdfsample[1:10]

## Compute atom pair library
ap <- sdf2ap(sdfset)

## ...
# fpset1024 <- names(rev(sort(table(unlist(as(apset, "list")))))[1:1024]) # uses 1024 most frequent APs
# apfpset <- desc2fp(x=apset, descnames=fpset1024, type="matrix")

}
\keyword{ utilities }
