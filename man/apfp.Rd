\name{apfp}
\alias{apfp}
\docType{data}
\title{
Frequent Atom Pairs
}
\description{
Ranked set of 4096 most frequent atom pairs observed in the compound collection from DrugBank with a MW < 1000. Their atom pairs were generated with the \code{sdf2ap} function. The provided data frame is sorted row-wise by atom pair frequency and only the 4096 most frequent atom pairs are included. This data set can be used as predefined atom pair selection when computing atom pair fingerprints with the \code{desc2fp} function.
}
\usage{data(apfp)}
\format{
  Object of class \code{data.frame}. First column contains atom pair (AP) IDs and the second column their frequency in DrugBank compounds.
}
\details{
Object stores 4096 most frequent atom pairs generated from DrugBank compounds.
}
\source{
DrugBank: http://www.drugbank.ca/
}
\references{
Chen X and Reynolds CH (2002). "Performance of similarity measures
in 2D fragment-based similarity searching: comparison of structural descriptors
and similarity coefficients", J Chem Inf Comput Sci.
}
\examples{
data(apfp)
apfp[1:4,]
}
\keyword{datasets}
