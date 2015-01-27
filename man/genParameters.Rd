\name{genParameters}
\alias{genParameters}
\title{
	Generate Parameters
}
\description{
	Generate statistics from a fingerprint database for use in caluclating z-scores, E-values, 
	and p-values later.
}
\usage{
	genParameters(fpset, similarity = fpSim, sampleFraction = 1, ...)
}
\arguments{
  \item{fpset}{
	  The database of fingerprints. Needs to be in the format expected by the similarity function.
	  For the default similarity function, this would be an \code{FPset}.
	}
  \item{similarity}{
	  A function to compute the similarity between two fingerprints. The first argument should be 
	  a single query and the second argument should be a set of fingerprints.
	}
  \item{sampleFraction}{
		The fraction of all pairs to use for estimating parameters. See Details section.
	}
  \item{\dots}{
		Extra parameters will be passed on to the similarity function.
	}
}
\details{
	A beta function will be fit to the distribution of similarity scores
	produced by the given similarity function. By default, all pairwise similarities will be
	computed. Since this can be expensive for large databases, one can also sample pairs to use.
	This can be done by setting \code{sampleFraction} to the fraction of all pairwise similarities to
	use. For example, for a database of 100 fingerprints, there are 10,000 pairs. Setting
	\code{sampleFraction} to 0.5 will result in only 5,000 pairs being used to estimate the
	parameters.

	Parameters are conditioned on the number of set bits. This function therefore groups
	fingerprints by the number of set bits they have and then estimates parameters for each group.
	A set of global parameters is also estimated and returned for use in cases where there was not
	enough data to estimate the parameters for a particular number of set bits.

}
\value{
	A data frame with  the following columns:
	\item{count}{ The number of similarities used to estimate these parameters}
	\item{avg}{ the mean}
	\item{variance}{ the variance}
	\item{alpha}{ The alpha paramber of the Beta function}
	\item{beta}{ The beta parameter of the Beta function}

	There will be a row for each possible count of 1 bits. So for a database of 1024 bit
	fingerprints, there will be 1025 rows for the possible values of 0-1024 bits. There will also be one
	additional row at the end with the global parameters. This can be used for cases where there
	are no parameters estimated for the current query 1-bit count.
}
\references{
	Pierre Baldi and Ramzi Nasr, "When is Chemical Similarity Significant? The Statistical
		Distribution of Chemical Similarity Scores and Its Extreme Values"
		Journal of Chemical Information and Modeling 2010 50 (7), 1205-1222
}
\author{
	Kevin Horan
}
\seealso{
	\code{\link{fpSim}}
}
\examples{

	library(ChemmineR)
	data(apset)
	fpset=desc2fp(apset) #get a fingerprint database

	params = genParameters(fpset)
	scores = fpSim(fpset[[1]],fpset,parameters=params,top=10)
}
\keyword{evalue}
\keyword{pvalue}
\keyword{fingerprint}
\keyword{similarity}
