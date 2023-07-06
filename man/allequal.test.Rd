\name{allequal.test}
\alias{allequal.test}
\title{all.equal with expected outcome test}
\description{
 Wrapper for all.equal, with ability to specify expected outcome.
}
\usage{
allequal.test(target,current,expect=TRUE)
}
\arguments{
  \item{target}{R object}
  \item{current}{R object}
  \item{expect}{logical, expected result of all.equal(x,y)}
}
\value{stop() if all.equal(x,y) != expect, else returns expect }

\author{Jonathan Wand}
\references{
  Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}

\examples{

a <- c(1:5)
b <- a+1

allequal.test(a,a,expect=TRUE)
allequal.test(a,b,expect=FALSE)
## these throw stop() but here we wrap in try() for example usage
try(allequal.test(a,a,expect=FALSE))
try(allequal.test(a,b,expect=TRUE))

}
\keyword{utilities}
