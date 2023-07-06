\name{anchors.chopit.fit}
\alias{anchors.chopit.fit}
\title{Compound Hierarchical Ordered Probit (CHOPIT)}
\description{
Fit function for pooled ordered probits
}
\usage{
  anchors.chopit.fit(data, parm, count, options )
}
\arguments{
    \item{data}{object of class anchors.data.  See help(anchors.data) for requirements.}
    \item{parm}{object of class anchors.chopit.parm.}
    \item{count}{object of class anchors.data.count.}
    \item{options}{object of class anchors.options.}
}
\details{

  This function is unlikely to be used directly by a user.
}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}
}

\seealso{\code{\link{anchors}},\code{\link{chopit}}}

\keyword{models}
