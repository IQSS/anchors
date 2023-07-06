\name{replace.list}
\alias{replace.list}
\title{Updating contents of one list using a second list}
\description{
 Update or insert named elements into a list using a second source list.
}
\usage{
replace.list( old, new )
}
\arguments{
  \item{old}{The target list to be updated}
  \item{new}{The soure list whose elements will be inserted into 'old'}
}

\value{An updated list will be returned -- note that the original }

\note{

    'old' is unchanged: remember to assign the result.

    Also, this is  recursive function if 'new' is a list of lists.
}

\author{Jonathan Wand}
\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}}}

\examples{

  ## replace y, and insert w
  a <- list( x = 1, y = 2, z = 3)
  b <- list( y = 4, w = 5)
  replace.list( a, b)

  ## recursive
  a <- list( x = list( y = 1, z = 2) )
  b <- list( x = list( y = 2, w = 3) )
  replace.list( a, b)

  ## if there is any disagreement between structure of old and new
  ## then structure of new list replaces structure of old list
  a <- list( x = 1, y = list( y=1,z=2))
  b <- list( x = list( y = 2, w = 3) , y = -9)
  replace.list( a, b)

}

\keyword{manip}
