\name{barplot.anchors.order}
\alias{barplot.anchors.order}
\title{Barplot of frequency of vignette orderings}
\description{Barplot of frequency of vignette orderings}
\usage{
  \method{barplot}{anchors.order}(height, ..., top=20)
}
\arguments{
  \item{height}{anchors.order object}
  \item{...}{arguments for \code{barplot}}
  \item{top}{Show 'top' most common combinations}
}
\details{

   Each vignette is represented by the integer indexing it's place as
   passed to the function.  For example, let

   fo <- list(vign = cbind(vign1,vign2,vign3))

   then 1=vign1, 2=vign2, 3=vign3. In the case of an anchors.order
   object produced with the option ties="interval", tied variables are
   represented by being included in brackets.

   z  <- anchors( fo, data, method="order",
                 options=anchors.options(ties="interval"))

  then a row with "\{1,3\},2" indicates that vign1 and vign3 are tied,
  and vign2 has a higher value than both of them.

}
\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}},\code{\link{anchors.order}}}
\examples{
## see help(anchors.order)  and demo(anchors.plot) for examples
}
\keyword{hplot}
