\name{anchors}
\alias{anchors}
\alias{anchors.rank}
\title{Non-parametric analysis of surveys with vignette anchors}
\description{
  Non-parametric analysis of surveys with vignette anchors.
}
\usage{
anchors(formula, data, method = c("B","C"),
        options=anchors.options(),
        subset, combn = FALSE, na.action = na.omit)
}
\arguments{
  \item{formula}{A list of named formulas giving a symbolic description
  of the model to be fit. See Details below.}
  \item{data}{A data frame or matrix.  See Details below.}
  \item{method}{Single string indicating method of analysis.  See
    Details below.}
  \item{options}{A list of class 'anchors.options', produced by function
    \code{anchors.options}; }
  \item{subset}{Logical expression indicating elements or rows to keep:
          missing values are taken as false; equivalent to function \code{subset}.}
  \item{combn}{Logical (TRUE/FALSE).  Calculate key statistics for all combinations of subsets of vignettes
    specified in formula=list(vign=cbind(...))}
  \item{na.action}{Specify the "NA action" which should be applied to
  'data'.  Default is R standard option 'na.omit'; can also be set to
  'na.fail'.}
}

\details{

     The formula syntax is designed to allow the same list of formulas
     to be used for both anchors and chopit.
     If a user would like to use exactly the same cases with
     \code{chopit} as with \code{anchors} (i.e., drop cases with any
     missing responses or missing \code{cpolr} values), and vise versa (i.e.,
     also drop any case with missing values in the \code{tau=},
     \code{tau1=} formulae), then the user should use the same list of
     formula for all methods and use the option
     \code{anchors.options(delete="maximal")}.  See example below.

      Only options relevant to \code{anchors} are discussed here.

      The named list() of formulas may include

      self: LHS self-response variable; RHS is simply 1 (needed to
      define a formula); actually anything can be put on the RHS, and it
      will not by default be used.

      vign: vignette responses, bound together by cbind (as in glm
    binomial syntax), with ONLY '~ 1'  as RHS of equation (formulas must
    have RHS).

    cpolr: OPTIONAL. linear predictors used by the censored ordered
    probit, \code{cpolr}, model to
    break ties in
    rank based models, method="B" or "C".  Default is cpolr = ~ 1 (i.e.,
    an intercept only).

    Example:

    \preformatted{
    fo <- list(self = xsayself ~ 1,
              vign = cbind(xsay3,xsay1)    ~ 1,
	      cpolr = ~ age )
	      }

    NOTE: one can also use a single formula as a short hand: specifying

    fo <- self ~ vign1 + vign2

    is equivalent to

    \preformatted{
    fo <- list(self = self ~ 1
               vign = cbind(vign1,vign2) ~ 1,
               cpolr= ~ 1)
	     }


  All of the response variables must be in the form of consequetive non-negative
  numeric integers, i.e., 1, 2, ... K.

    ** \code{anchors} currently does not support factor responses.**

  The method="B" and "C" requires that cases with any response
  that is missing be dropped.  Any cases with missing values
  the covariates specified in the \code{cpolr=} formula are also
  dropped.  BUT by default these methods ignore missing values in the
  covariates specified by \code{tau=}, \code{tau1=}, or the right hand
  side of \code{self=}.
}

\value{

  Specific values listed below are meant for programmers; data analysts
  should see the values produced by the \code{summary.anchors} function.

  \item{rank}{List object containing:

    \describe{
      \item{span}{ A matrix with columns  rows
	corresponding to the observations in the data for which both the
	self-assessment and the vignettes are observed.

	The first column contains the lower bound of the rank interval, and
	the second column contains the upper bound.  If a row has the same
	value in both columns, then the rank value is a scalar value,
	otherwise the rank lies somewhere in an interval.

	The columns will be named \code{Cs} and \code{Ce} for type "C",
	\code{Bs} and \code{Be} for type "B".
      }

      \item{weight}{matrix of weights, inverse proportion to span of interval
	values}
      \item{max}{Maximum rank possible for non-parametric method

	For "C" max is 2J+1, where J is the number of vignettes.

	For "B" max is J+1.
      }

      \item{n.interval}{Number of cases that have interval values of C}
  }} %% end of rank

  \item{summary}{ list of summary statistics about ranks}

  \item{minentropy}{ Matrix specifying how interval rank values would be
    allocated to scalar rank values to minimize entropy.}

  \item{cpolr}{ a \code{cpolr} model that was used to estimate how to
    allocate interval rank values to scalar rank values}

  \item{type}{ The type of nonparamtric analysis (B or C) }

}




\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

Wand, J. (2013), Credible Comparisons Using Interpersonally Incomparable Data: Nonparametric Scales with Anchoring Vignettes. \emph{American Journal of Political Science}, 57: 249-262. \doi{10.1111/j.1540-5907.2012.00597.x}

King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 15(1), 46-66. \doi{10.1093/pan/mpl011}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}

}
\seealso{
  \code{\link{anchors.order}},
  \code{\link{anchors.options}},
  \code{\link{chopit}},
}

\examples{
data(mexchn)

fo <- list(self = xsayself ~ age,
           vign = cbind(xsay3,xsay1)    ~ 1,
           tau  = ~ age + china,
           tau1 = ~ age + china + male,
           cpolr= ~ china + male
    )

a2 <- anchors(fo, data = mexchn, method="B")
summary(a2)

a3 <- anchors(fo, data = mexchn, method="C")
summary(a3)


}
\keyword{models}
