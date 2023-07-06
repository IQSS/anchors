\name{fitted.anchors.cpolr}
\alias{fitted.anchors.cpolr}
\title{Conditional and unconditional prediction for censored ordered probit}
\description{
  Conditional and unconditional prediction for censored ordered
  probit.  Unconditional prediction returns the fitted values (predicted
probabilities) from the \code{\link{cpolr}} object.  Conditional prediction
takes the observed range of the diff-corrected self-response output from
\code{\link{anchors}} and renormalizes the predicted
probabilities for each observation.}

\usage{
   \method{fitted}{anchors.cpolr}(object, average = FALSE, unconditional = FALSE, ...)
}
\arguments{
  \item{object}{anchors.cpolr object}
  \item{average}{a logical value.  See \code{values} below for more details.}
  \item{unconditional}{Set to TRUE if you submit an
    anchors.object AND want the unconditional probabilities returned.
    One case that you would submit a anchors.rank object is
  if you did subsetting for the anchors object but not for the cpolr
  object, and want the intersection of the two objects used for the
  unconditional probabilities.}
  \item{\dots}{required for S3, but any other options will be ignored.}
}
\value{If \code{average = FALSE}, a matrix of predicted probabilities
  with rows corresponding to observations, and columns corresponding to
  categories.

  If \code{average = TRUE}, the matrix of predicted probabilities
  (conditional or unconditional) is summarized to a vector (summed by categories,
  then renormalized to sum to 1).

  If \code{anchors} object has been specified, then each observation is
  renormalized to fall into the range of the diff-corrected
  self-response for that observation.  If there are no ties for a given
  observation, then that observation is a
  vector consisting of (k-1) zeros and 1 one.  If there are ties, then
  the predicted probabilities for that observation are renormalized to
  fall within the diff-corrected range.

  If \code{anchors} object is omitted, identical to the matrix of predicted
  probabilities from the \code{\link{cpolr}} output.

}

\author{Jonathan Wand}

\references{
Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S.}  4th edition.  Springer.

Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}

}
\seealso{\code{\link{anchors}}, \code{\link{cpolr}}}
\examples{

## see examples in anchors

}

\keyword{models}
