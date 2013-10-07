\name{fingerprintOB}
\alias{fingerprintOB}
\title{
	Fingerprints from OpenBabel
}
\description{
	Generates fingerprints from SDFsets using OpenBabel. The name of the fingerprint can
	also be set and can be anything available through OpenBabel. You can see what this list is by
	executing "obabel -L fingerprints".  Results are returned as an FPset.
}
\usage{
fingerprintOB(sdfSet, fingerprintName)
}
\arguments{
  \item{sdfSet}{
	  Input compounds to generate fingerprints for.
	}
  \item{fingerprintName}{
	  The name of the fingerprint in Open Babel. A list of available names can be found
	  by executing "obabel -L fingerprints".  Currently that list is: "FP2", "FP3", "FP4", and "MACCS". 
	}
}
\value{
	An FPset with an element for each given compound.
}
\author{
	Kevin Horan
}
\examples{
	\dontrun{
		data(sdfsample)
		fpset = fingerprintOB(sdfsample)
	}
}
