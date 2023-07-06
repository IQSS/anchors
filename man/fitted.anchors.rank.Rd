\name{fitted.anchors.rank}
\alias{fitted.anchors.rank}
\title{Fitted values of non-parametric models}
\description{Fitted values of non-parametric models}
\usage{
   \method{fitted}{anchors.rank}(object, ..., ties = c("omit", "uniform", "cpolr",
                 "minentropy"), average = FALSE, unconditional = FALSE)
}
\arguments{
  \item{object}{object of class \code{\link{anchors.rank}}.}
  \item{\dots}{required for S3, but any other options will be ignored.}
  \item{ties}{method for allocating ties/interval valued cases to a
    scalar value}
  \item{average}{a logical value.  See \code{values} below for more details.}
  \item{unconditional}{Set to TRUE if you submit an
    anchors.object AND want the unconditional probabilities returned.
    One case that you would submit a anchors.rank object is
  if you did subsetting for the anchors object but not for the cpolr
  object, and want the intersection of the two objects used for the
  unconditional probabilities.}
}
\value{If \code{average = FALSE}, a matrix of predicted probabilities
  with rows corresponding to observations, and columns corresponding to
  categories.

  If \code{average = TRUE}, the matrix of predicted probabilities
  (conditional or unconditional) is summarized to a vector.

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

## Basic usage: see anchors

}

\keyword{models}
