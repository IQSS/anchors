\name{plot.anchors.combn}
\alias{plot.anchors.combn}
\title{Plot results from anchors( method='entropy') }
\description{Plot results from anchors( method='entropy') }
\usage{
  \method{plot}{anchors.combn}(x, ... , xy )
}
\arguments{
   \item{x}{anchors.combn object}
   \item{...}{additional options for plot}
   \item{xy}{name two columns from anchors.combn to plot (case
     insensitive, and you only need to specify first three letters of
     each name); must be a      vector of strings of length 2.}
}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

}
\seealso{\code{\link{anchors}},\code{\link{anchors.combn}}}
\examples{
## see help(anchors.combn) and demo(anchors.plot) for examples
}
\keyword{hplot}
