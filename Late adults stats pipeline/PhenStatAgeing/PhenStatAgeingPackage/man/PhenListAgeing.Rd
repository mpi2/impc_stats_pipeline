\name{PhenListAgeing}
\alias{PhenListAgeing}
\title{
	Method "PhenListAgeing"
}
\description{
	Function to create PhenListAgeing object from a PhenList object.
}
\usage{
	PhenListAgeing(
		PhenListobject,
		DOE = NULL,
		DOB = NULL,
		d.threshold = 16 * 7,
		debug = TRUE
	)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PhenListobject}{
		Instance of the PhenList class; mandatory argument
}
  \item{DOE}{
  The name of the column for the date_of_experiment in the PhenList object. If left NULL then the input PhenlistObject will be returned. Default NULL
}
  \item{DOB}{
  The name of the column for the date_of_birth in the PhenList object. If left NULL then the input PhenlistObject will be returned. Default NULL
}
  \item{d.threshold}{
  The threshold in age (DOE-DOB)  to specify LifeStage early/late levels. The function uses as.Date(DOE)-as.Date(DOB) to calculate the age of the samples. The default is 16 weeks (16*7 days)
}
  \item{debug}{
  A logical flag. Set to TRUE to see debug messages. Default TRUE
}
}

\value{
	Provided DOE and DOB are not NULL,  a PhenListAgeing object that is similar to PhenList object with an extra column called LifeStage with two levels Early and Late. Otherwise,  the output is similar to the input PhenList object.
}

\author{
Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{testDatasetAgeing}}
}
\examples{
library(PhenStat)
file = system.file("extdata", "test1.csv", package = "PhenStat")

####################################################################
# PhenListObject
####################################################################
PhenListObject = PhenList(dataset = read.csv(file, na.strings = '-'),
 testGenotype = "Sparc/Sparc")
class(PhenListObject)

####################################################################
# PhenListAgeing object
####################################################################
PhenListAgeing   = PhenListAgeing(
	PhenListobject = PhenListObject,
	DOE = 'Batch',
	DOB = 'Birth.Date',
	d.threshold = 99
)
plot(PhenListAgeing)
class(PhenListAgeing)
				
}

\keyword{ PhenListAgeing }

