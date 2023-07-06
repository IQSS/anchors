\name{trim.data}
\alias{trim.data}
\title{Trim a dataset to have same cases present in an anchors.data object}
\description{
  Given a dataset X that was used in the process of creating an anchors.data
  object Y (which may have had cases deleted due to missing values),
  trim cases of X to include only cases that would be used in an
  anchors() analyis of Y.  Matching is done by rownames.
}
\usage{
trim.data( data, anchors )
}
\arguments{
  \item{data}{A data.frame or matrix}
  \item{anchors}{An object returned by \code{anchors} }
}

\value{A data.frame or matrix with cases that match rows of
  anchors$data$z0; matching done by rownames }

\note{

    'data' is unchanged: remember to assign the result.

}

\author{Jonathan Wand}
\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}}}

% \examples{
%
% }

\keyword{manip}
