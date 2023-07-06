\name{replace.value}
\alias{replace.value}
\title{Replaces occurences of a value with another value in set of columns}
\description{
 Replaces a single value in a set of columns
 with another given value.   This makes it easy to change the default
 missing value indicator, for example.
}
\usage{
replace.value( data, names, from=NA, to=as.integer(0), verbose = FALSE)
}
\arguments{
  \item{data}{data frame}
  \item{names}{a vector of character strings identifying columns to be updated }
  \item{from}{value to find and replace, can specify vectors}
  \item{to}{replacement value, must be a scalar}
  \item{verbose}{prints warnings if typeof(from) != typeof(to)}
}

\value{ A data frame, with the same number of observations as the input
  \code{data}, but with replaced values as specified.}

\details{
  Will replace storage mode of 'to' with mode of 'from' in dataset.
}

\author{Jonathan Wand}
\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}}}
\examples{
## data with zeroes as missing values in responses
data(poleff)
## data with NA missing values in responses
data(poleffna)

## convert NA to 0:
dd  <- replace.value(poleffna,c("xsayself","xsay1","xsay2","xsay3","xsay4","xsay5"))

## convert 0 to NA:
dd2 <- replace.value(poleff,c("xsayself","xsay1","xsay2","xsay3","xsay4","xsay5"),0,as.double(NA))

}
\keyword{manip}
