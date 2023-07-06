\name{fitted.cpolr}
\alias{fitted.cpolr}
\title{Conditional and unconditional prediction for censored ordered probit}
\description{
  Conditional and unconditional prediction for censored ordered
  probit.  Unconditional prediction returns the fitted values (predicted
probabilities) from the \code{\link{cpolr}} object.  Conditional prediction
takes the observed range of the diff-corrected self-response output from
\code{\link{anchors}} and renormalizes the predicted
probabilities for each observation.}

\usage{
   \method{fitted}{cpolr}(object, anchors, average = FALSE, unconditional = FALSE, ...)
}
\arguments{
  \item{object}{output from \code{\link{cpolr}}.}
  \item{anchors}{leave missing for unconditional prediction (or set unconditional=TRUE).  For
    conditional prediction, specify the object of class anchors.rank
    used to run cpolr originally.}
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
\note{If the \code{anchors} object is made using a subset of the data
  used to create the \code{cpolr} object, then invoking
  \code{fitted.cpolr} will not use the same cases in calculating the
  conditional probabilities as it would if the \code{anchors} object is
  omitted!

  If you want to have the same cases used in the unconditional
  calculation as in the conditional with a subsetted \code{anchors} object,
  then include  \code{anchors} object and set
  \code{unconditional.override = TRUE}.

  }

\author{Jonathan Wand}

\references{
  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
     Statistics with S.}  4th edition.  Springer.

Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 15(1), 46-66. \doi{10.1093/pan/mpl011}

}
\seealso{\code{\link{anchors}}, \code{\link{cpolr}}}
\examples{

## Basic usage: see cpolr

}

\keyword{models}
