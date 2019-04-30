\name{plot.SmoothWin}
\alias{plot.SmoothWin}

\title{
 Plot function for the SmoothWin object
}
\description{
  This function plots a SmoothWin object
}
\usage{
  \method{plot}{SmoothWin}(x, ylab = 'Response', col= NULL, digits = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  SmoothWin object
}
  \item{ylab}{
  Label on the y axis
}
  \item{col}{
   Colour parameter for the points. Set to NULL to use the default colouring (spectrum colouring)
  }
  \item{digits}{
   The number of visible digits for l, k and SWS 
  }
  \item{\dots}{
  Optional parameters that can be passed to the `plot`/`qqPlot` function. See `car` package for the qqPlot function
}
}
\author{
	Hamed Haselimashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
	\code{\link{SmoothWin}}
}

\examples{
	example(SmoothWin)
}
