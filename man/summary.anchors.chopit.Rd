\name{summary.anchors.chopit}
\alias{summary.anchors.chopit}
\title{Summary of CHOPIT Analysis}
\description{
  Prints estimates and -log-likelihood from the out of CHOPIT model.
}
\usage{
 \method{summary}{anchors.chopit}( object, ..., digits = 4 )
}
\arguments{
  \item{object}{ The object of class "anchors.chopit", as returned by
    \code{\link{anchors}}}
  \item{...}{additional arguments affecting the summary produced}
  \item{digits}{Number of digits you would like to have printed, using
    \code{\link{round}} function.}
}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}
}

\seealso{\code{\link{anchors}},\code{\link{chopit}}}
\examples{
## see help(anchors)
}
\keyword{models}
