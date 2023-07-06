\name{anchors.data}
\alias{anchors.data}
\title{Organized data from surveys with anchoring vignettes}
\description{
  Create data object that will be used by anchors().
}
\usage{
anchors.data(formula,
             data,
             method,
             subset,
             na.action   = na.omit,
             na.response = c(NA, 0),
             min.response = 1,
             delete = c("minimal","maximal"), debug=0)
}
\arguments{
  \item{formula}{A list of named formulas giving a symbolic description
    of the model to be fit. See help(anchors) for details.}

  \item{data}{A data frame or matrix}

  \item{method}{Single string indicating method of analysis.
    See help(anchors) for details.}

  \item{subset}{Logical statement as used by function \code{subset}.}

  \item{na.action}{a function which indicates what should happen when the data
    contain 'NA's.  Note also details of 'delete' in \code{anchors.options}.}

  \item{na.response}{vector of numeric values that should also be
    considered missing in self or vignette responses.  Default is c(0,NA). }

  \item{min.response}{Default=1.  This is a sanity check value.  The
    code is designed to make sure that the user has not passed in values
    that are less than 1 and not included in na.response list.}

  \item{delete}{"minimal" deletes only cases with missing values that
    affect component of model

    "maximal" forces listwise deletion on the basis
    of ALL variables in the entire formula list EVEN if vars not used by
    method
  }

  \item{debug}{Default: 0}

}

\value{
  Return function will be of class 'anchors.data'.
}

\details{

  All of the response variables must be in the form of consequetive non-negative
  numeric integers, i.e., 1, 2, ... K.

  ** \code{anchors} currently does not support factor responses.**

  The method="chopit" model does *not* require that every
  respondent answer all vignette or the self questions to be included in
  the model; only answers that are non-missing are used in the
  likelihood function.

  The method="B", "C", and "order" requires that cases with any response
  that is missing be dropped.

  By default, \code{anchors} only deletes those cases with missing
  values that affect the method of analysis requested.  For example,

  If a user would like to use exactly the same cases with
  method="chopit" as "B", "C" and "order" (i.e., drop cases with any
  missing responses), and vise versa (i.e., also drop any case with
  missing values in the \code{tau=}, \code{tau1=}, or \code{self=}
  formula), then the user should use the same list of formula for all
  methods and use the option \code{anchors.options(delete="maximal")}.

}
\section{anchors() and anchors.options()}{

  Specifying no options is equivalent to

  \code{anchors(..., anchors.options(delete = "minimal")) }

  Alternative values include,

    \describe{
      \item{\code{delete}}{Choose one from:

	"minimal": delete only cases with missing values in components
	needed for current anchors analysis

	"maximal": forces listwise deletion on the basis
        of any missing value in ALL variables in the formula list EVEN
        if not vars not used by method.  Also deletes cases with ANY
	missing values in any responses (self or vignettes).  This
	enables method="chopit" to use the same data as used by
	method="B" or "C", and vice versa.


      }
    }

}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}} }
\examples{
## load data and make a copy
data(mexchn)
dta <- mexchn

## insert missing values into covariates for the purpose of this demonstration
dta <- replace.value( dta, "educyrs", from=1:9  , to = NA)
dta <- replace.value( dta, "age"    , from=30:40, to = NA)

## formula that will be used throughout
fo <- list(self = xsayself ~ male + educyrs,
           vign = cbind(xsay3, xsay1) ~ 1,
           tau  = ~ educyrs)

## 'C' uses only cases with no missing responses (self, vign)
a1 <- anchors( fo, dta, method="C")
## number of cases used:
##  y0 = number of cases with self-responses
##  z0 = number of cases with all vignette responses
unlist(lapply(a1$data, NROW))[1:2]

## 'chopit' keeps case with some missing responses
##  but drops others with missing covariates
a2 <- chopit( fo, dta)
## number of cases used:
##  y0 = number of cases with self-responses
##       AND no missing covariates
##  z0 = number of cases with AT LEAST ONE vignette response
##       and no missing covariates
unlist(lapply(a2$data, NROW))[1:2]

## with ' delete = "maximal" '
## both procedures use the exact same cases
## (IF the same formula is used for both methods)

## y0 and z0 here
a3 <- anchors( fo, dta, method="C", anchors.options(delete="maximal"))
unlist(lapply(a3$data, NROW))[1:2]

## is the same as y0 and z0 for chopit:
a4 <- chopit( fo, dta, options=anchors.options(delete="maximal"))
unlist(lapply(a4$data, NROW))[1:2]

}

\keyword{models}
