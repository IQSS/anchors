\name{anchors.order}
\alias{anchors.order}
\alias{anchors.order.calc}
\title{Calculate frequency of vignette orderings}
\description{Calculate frequency of vignette orderings}
\usage{
  anchors.order(formula, data, ties = c("set", "nominal", "random",
                 "mset"), subset, na.action = na.omit)
}
\arguments{
  \item{formula}{A list of named formulas giving a symbolic description
    of the model to be fit. See Details below.}

  \item{data}{A data frame or matrix.  See Details below.}

  \item{ties}{"set": Groups ties as sets

    "nominal":  Breaks ties by order of vignettes given

    "random": Breaks ties randomly }

  \item{subset}{Logical expression indicating elements or rows to keep:
    missing values are taken as false; equivalent to function \code{subset}.}

  \item{na.action}{Specify the "NA action" which should be applied to
    'data'.  Default is R standard option 'na.omit'; can also be set to
    'na.fail'.}

}

\details{
  The formula syntax is designed to allow the same list of formulas
     to be used for both anchors and chopit and anchors.order.
     If a user would like to use exactly the same cases with
     \code{chopit} as with \code{anchors.order} (i.e., drop cases with any
     missing responses or missing \code{cpolr} values), and vise versa (i.e.,
     also drop any case with missing values in the \code{tau=},
     \code{tau1=} formulae), then the user should use the same list of
     formula for all methods and use the option
     \code{anchors.options(delete="maximal")}.  See example below.

      Only options relevant to \code{anchors.order} are discussed here.

      The named list() of formulas may include

      self: LHS self-response variable; RHS is simply 1 (needed to
      define a formula); actually anything can be put on the RHS, and it
      will not by default be used.

      vign: vignette responses, bound together by cbind (as in glm
    binomial syntax), with ONLY '~ 1'  as RHS of equation (formulas must
    have RHS).

    Example:

    \preformatted{
    fo <- list(self = xsayself ~ 1,
              vign = cbind(xsay3,xsay1)    ~ 1)
	      }

    One can also use a single formula as a short hand: specifying

    fo <- self ~ vign1 + vign2

    is equivalent to

    \preformatted{
    fo <- list(self = self ~ 1
               vign = cbind(vign1,vign2) ~ 1)
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

   Each vignette is represented by the integer indexing it's place as
   passed to the function.  For example, let

   fo <- list(vign = cbind(vign1,vign2,vign3))

   then 1=vign1, 2=vign2, 3=vign3. In the case of an anchors.order
   object produced with the option ties="set", tied variables are
   represented by being included in brackets.

   \preformatted{
   anchors.order( fo, data, ties="set")
   }

  then a row with "\{1,3\},2" indicates that vign1 and vign3 are tied,
  and vign2 has a higher value than both of them.


  \item{labels}{Description of order of vignettes}
  \item{freq}{Frequency of order of vignettes}
}

\author{Jonathan Wand and Dan Hopkins}
\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 15(1), 46-66. \doi{10.1093/pan/mpl011}
}
\seealso{\code{\link{anchors}},\code{\link{barplot.anchors.order}}}
\examples{
data(mexchn)
z <- anchors.order( ~ xsay1 + xsay2 + xsay3 + xsay4 + xsay5, mexchn, ties  = "set")
summary(z,top=10,digits=3)
barplot(z)

z <- anchors.order( ~ xsay5 + xsay4 + xsay3 + xsay2 + xsay1, mexchn, ties  = "set")
summary(z,top=10,digits=3)

## other option
z <- anchors.order( ~ xsay5 + xsay4 + xsay3 + xsay2 + xsay1, mexchn, ties  = "nominal")
summary(z,top=10,digits=3)
barplot(z)


}
\keyword{hplot}
