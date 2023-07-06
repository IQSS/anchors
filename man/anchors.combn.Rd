\name{anchors.combn}
\alias{anchors.combn}
\title{Calculate known minimum or estimated entropy for survey vignettes}
\description{
  Calculate known minimum or estimated entropy for survey vignettes
}
\usage{
   anchors.combn(adata, fdata, type, options)
}
\arguments{
  \item{adata}{object of class anchors.data}
  \item{fdata}{data frame from which adata was built }
  \item{type}{either "B" or "C", specifying nonparametric method}
  \item{options}{object of class anchors.options}
}
\value{A list of class \code{anchors.combn}, containing elements
 \item{vign.rank}{a matrix with columns corresponding to the vignette
    selection, the estimated entropy (if covar is specified), the
  known entropy, the number of cases with intervals for that particular vignette
  selection, and 2J+1 number of categories for the given
  vignettes.}
 \item{N}{the number of observations after list-wise deletion}
 \item{self}{a character string containing the variable name of the
    self-response question}
 \item{vign}{a character vector containing the variable names of the
      vignette questions.  The first element of the vector corresponds
      to vignette 1, the second vignette 2, etc. }
  }

\author{Jonathan Wand with Dan Hopkins and Olivia Lau}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}
   }

\seealso{\code{\link{anchors}}}


\examples{

  ## ra <- anchors( fo, data, method="C", combn=TRUE)
  ## summary( ra , combn = "Minimum" )
  ## plot( ra, type="xy" , xy =  c("minimum","intervals") )

}

\keyword{htest}
