\name{summary.anchors.order}
\alias{summary.anchors.order}
\title{Calculate frequency of vignette orderings}
\description{Calculate frequency of vignette orderings}
\usage{
  \method{summary}{anchors.order}(object, top, digits = getOption("digits"),verbose=FALSE, ...)
}
\arguments{
   \item{object}{anchors.order object}
   \item{top}{Show 'top' most common combinations}
   \item{digits}{digits to print}
   \item{verbose}{additional summaries}
   \item{...}{additional arguments}
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
## see example in help(anchors.order)
}
\keyword{hplot}
