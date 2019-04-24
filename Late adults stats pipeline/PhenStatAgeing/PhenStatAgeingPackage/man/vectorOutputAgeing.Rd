\name{vectorOutputAgeing}
\alias{	vectorOutputAgeing}

\title{
	Method "vectorOutputAgeing"
}
\description{
	Wrapper for the model fitting results. Returns model fitting and testing results in a list or JSON format.
}
\usage{
	vectorOutputAgeing(
		object,
		othercolumns = NULL,
		JSON  = FALSE,
		Null  = FALSE,
		...
	)
}
%-maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
		Instance of the testDatasetAgeing result;	mandatory argument
	}
	\item{othercolumns}{
		A list of column names that must be included in the results. Default NULL
	}
	\item{JSON}{
		A logical flag. Setting to TRUE for the JSON output otherwise, the function returns a list
	}
	\item{Null}{
		A logical flag. Setting to TRUE forces the function to return results even if the testDatasetAgeing returns a failure message
	}
	\item{\dots}{
		Other parameters that can be passed to toJSON() function in "jsonlite" library
	}
}

\value{
	A list of values or a JSON object depends on the "JSON" parameter
}

\author{
	Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{testDatasetAgeing}}
}
\examples{
	example(testDatasetAgeing)
}

