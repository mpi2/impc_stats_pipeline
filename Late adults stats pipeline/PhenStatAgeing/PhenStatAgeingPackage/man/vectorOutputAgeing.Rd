\name{vectorOutputAgeing}
\alias{	vectorOutputAgeing}

\title{
	Method "vectorOutputAgeing"
}
\description{
	Wrapper for the output of 'testDatasetAgeing'. Returns model fitting and results in a list or JSON format (StatPacket).
}
\usage{
vectorOutputAgeing(
	object                   ,
	othercolumns = NULL      ,
	JSON = FALSE             ,
	RemoveNullKeys = FALSE   ,
	ReportNullSchema = FALSE ,
	...
)
}
%-maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{
		'mandatory argument'. An instance of the testDatasetAgeing result
	}
	\item{othercolumns}{
		A list of column names that must be included in the results. Default NULL
	}
	\item{JSON}{
		Logical flag. Setting to TRUE for the JSON (StatPacket) output otherwise, the function returns a list
	}
	\item{RemoveNullKeys}{
		Logical flag. Setting to TRUE will remove all NULL elements from the output. Default is FALSE
	}
	\item{ReportNullSchema}{
		logical flag. Setting to TRUE forces the function to return results even if the testDatasetAgeing returns a failure message
	}
	\item{\dots}{
		Other parameters that can be passed to `toJSON()` function in the "jsonlite" library
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

