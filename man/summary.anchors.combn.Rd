\name{summary.anchors.combn}
\alias{summary.anchors.combn}
\title{Summary of anchors(..., method='entropy') analysis}
\description{Summary of anchors(..., method='entropy') analysis}
\usage{
 \method{summary}{anchors.combn}( object, ...,
                   sort = c("max","estimated","minimum","interval","span"),
                                    digits=3)
}
\arguments{
  \item{object}{ The object of class "anchors.combn", as returned by
    \code{\link{anchors}}}
  \item{...}{additional arguments affecting the summary produced.}
  \item{sort}{column by which to sort results}
  \item{digits}{sig digits to print}
}
\author{Jonathan Wand}
\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}

King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 15(1), 46-66. \doi{10.1093/pan/mpl011}

}
\seealso{\code{\link{anchors}}}
\examples{
## see anchors
}
\keyword{models}
