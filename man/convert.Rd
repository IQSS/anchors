\name{convert}
\alias{convert}
\title{Convert factor or character variables into integers}
\description{
  Use \code{convert} to transform factor or character variables into
  ordered integer values prior to using \code{anchors} or
  \code{chopit}.
}
\usage{
convert(vars, data, order, \dots)
}
\arguments{
  \item{vars}{a character vector identifying the variables in
    \code{data} to be converted to integers.}
  \item{data}{data frame containing \code{vars}.}
  \item{order}{a character vector identifying the factor levels to be
    converted to integers, in increasing order.  If \code{NA} is
    omitted here (default), \code{NA} values are retained as \code{NA}.
    If \code{NA} is included, then \code{NA} values are converted to the
    the last category (which can be adjusted via \code{\dots}.}
  \item{\dots}{addtional arguments passed to \code{factor}.}
}

\value{ A data frame, with the character or factor variables identified in \code{vars}
  replaced by integers.}

\author{Olivia Lau}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{
  \code{\link{anchors}},  \code{\link{chopit}}
}


\examples{

set.seed(391)
a <- sample(c("net gain", "little effect",  "modest loss", "major loss"),
            20, replace=TRUE)

dta <- as.data.frame( list(ncost=a) , stringsAsFactors  = FALSE)
print(dta)

newdat <- convert("ncost", data = dta,
                  order = c("net gain", "little effect",  "modest loss", "major loss"))

print(newdat)

}
\keyword{manip}
