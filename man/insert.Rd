\name{insert}
\alias{insert}
\title{Insert DIF-corrected variable into original data frame}
\description{
 Insert the DIF corrected variable into the original data frame, with
 missing values for observations for which it was impossible to
 calculate DIF correction (due to missingness in either the
 self-response or one or more of the vignette responses).
}
\usage{
insert(data, obj, vnames, overwrite, debug)
}
\arguments{
  \item{data}{Original data frame}

  \item{obj}{anchors.rank class object}

  \item{vnames}{A list of character vectors of length 2,which represent the
    variable names (in the output data frame) for the lower and upper
    bounds, respectively, of  the DIF corrected variables.

    By default, this will be derived from column names of B and C in
    the anchors.rank object: \code{B=c("Bs", "Be"),C=c("Cs", "Ce")}. }

  \item{overwrite}{A logical value (defaults to \code{FALSE}).
    Indicates whether to overwrite variables with names the same as
    \code{colnames} already exist in \code{data}. }

  \item{debug }{Default: 0}
}

\value{ A data frame, with the same number of observations as the input
  \code{data}, but with two or four additional variables, corresponding
  to the lower and upper bounds of the  DIF corrected variable.}

\author{Olivia Lau}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}},\code{\link{cpolr}} }


\examples{

data(freedom)
ra <- anchors(self ~ vign1 + vign3 + vign6, data = freedom, method="B")
freedom <- insert(freedom, ra )
names(freedom)

}
\keyword{manip}
