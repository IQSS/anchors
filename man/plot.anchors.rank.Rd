\name{plot.anchors.rank}
\alias{plot.anchors.rank}
\title{Plot components of anchors.rank objects}
\description{Plot components of anchors.rank objects.  Currently this
  only makes a scatter plot if anchors(...,combn=TRUE).  For plots of
  densities of ranks, see \code{barplot.anchors.rank} }
\usage{
  \method{plot}{anchors.rank}(x, ..., xy)
}
\arguments{
   \item{x}{anchors.rank object}
   \item{...}{for plot}
   \item{xy}{name two columns from anchors.combn to plot (case
     insensitive, and you only need to specify first three letters of
     each name); must be a vector of strings of length 2; see Details below.}
 }
 \details{

   Each observation in the plot is a subset of vignettes, represented by
   their index values as originally passed by the formula (see
   \code{anchors.order} for an example/description).

   To define the axis of the scatter plot, the columns to choose from
   include:

   "estimated": estimated entropy from cplor model

   "minimum": minimum entropy from minimum.entropy calculation

   "interval": number of cases with interval (non-scalar) ranks

   "span" average span of ranks (including scalar cases)

   "max" maximum B or C value for the subset of vignettes

   Choose two of the above, e.g., xy=c("minimum","interval")

   }


\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{
  \code{\link{anchors}},
  \code{\link{anchors.order}}
}
\examples{
## see help(anchors) and demo(anchors.plot) for examples
}
\keyword{hplot}
