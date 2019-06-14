\name{plot.PhenListAgeing}
\alias{plot.PhenListAgeing}

\title{
	plot for a PhenListAgeing object
}
\description{
	This function visualises a PhenListAgeing object
}
\usage{
	\method{plot}{PhenListAgeing}(x, vars = NULL, ...)
}

\arguments{
	\item{x}{
		PhenListAgeing object 
	}
	\item{vars}{
		Variable(s) of interest. The default is 'Batch', 'Genotype', 'Sex and 'LifeStage' if exists in the data
	}
	\item{\dots}{
		Optional parameters that can be passed to `Hmisc::plot.describe()`
	}
}


\author{
	Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{PhenListAgeing}}
}

\examples{
	example(PhenListAgeing)
}
