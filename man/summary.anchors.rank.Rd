\name{summary.anchors.rank}
\alias{summary.anchors.rank}
\title{Summary of non-parameteric anchors analysis}
\description{
  Summary of non-parameteric analysis of a set of vignettes and a self-assessment.
}
\usage{
 \method{summary}{anchors.rank}( object, ..., ties = c("omit", "uniform", "cpolr",
                 "minentropy"), combn = TRUE, digits = 3)
}
\arguments{
  \item{object}{ The object of class "anchors.rank", as returned by
    \code{\link{anchors}}}
  \item{...}{additional arguments affecting the summary produced.}
  \item{ties}{(optional) shows distribution of ranks by method for allocating ties/interval valued
    cases to a scalar value; can be one or more of the following:

     "uniform": allocates in equal proportions among scalar values

     "omit": drops cases with interval values

     "cpolr": uses censored ordered probit model to allocate ties

     "minentropy": allocates ties to produce minimum entropy

  }
  \item{combn}{print summary of combn if avaiable in anchors.rank object}
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
