\name{PhenListAgeing}
\alias{PhenListAgeing}
\title{
	Converting PhenList object into PhenListAgeing
}
\description{
	This function converts a PhenList object from PhenStat package into a PhenListAgeing object
}
\usage{
	PhenListAgeing(object, DOE = "Batch", DOB = "date_of_birth", d.threshold = 16 * 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
		PhenList object (for more details about the PhenList object see the manual in PhenStat package)
	}
	\item{DOE}{
		Date of experiment index. The default is 'Batch'
	}
	\item{DOB}{
		Date of birth index. The Default is 'date_of_birth'
}
  \item{d.threshold}{
  The threshold to discriminate young and old in days. The default is 16 weeks (16*7 days)
}
}


\author{
Hamed Haseli <hamedhm@ebi.ac.uk>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{testDatasetAgeing}}
}
\examples{
	example(testDatasetAgeing)
}

