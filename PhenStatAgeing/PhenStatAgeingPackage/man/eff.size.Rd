\name{eff.size}
\alias{eff.size}

\title{
	Effect size from lme or gls object
}
\description{
	This function computes the Cohen effect size from an lme or gls object
}
\usage{
	eff.size(object, data, depVariable, effOfInd = "Genotype")
}

\arguments{
	\item{object}{
		lme or gls object
	}
		\item{data}{
		a data frame containing the variables in the model. If leave Null then the function extracts the data from the model
	}
	\item{depVariable}{
		dependent variable 
	}
	\item{effOfInd}{
		effect of the index, that is the name of the variable. The default is 'Genotype'
	}
}
\details{
	The effect size is computed by \deqn{(\mu_1-\mu_2)/\sigma_{residuals}} where \deqn{\mu_1} and \deqn{\mu_2} are means of the response split on the levels of effOfInd.
}
\value{
	The effect size
}

\author{
	Hamed Haseli <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{testDatasetAgeing}}
}
\examples{
	example(testDatasetAgeing)
}
