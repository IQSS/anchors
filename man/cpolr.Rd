\name{cpolr}
\alias{cpolr}
\title{Censored ordered probit}
\description{
  Censored ordered probit for analysis of anchoring vignettes.  Used
  in the context of anchoring vignettes as a parametric model for
  breaking ties/interval in non-parametric ranks.
}
\usage{
  cpolr(formula, data, weights, start, ..., subset, na.action,
                 contrasts = NULL, Hess = TRUE, model = TRUE, method =
                 c("probit", "logistic", "cloglog", "cauchit"), debug = 0)
}

\arguments{
  \item{formula}{A formula representing 'C' range produced by
    \code{anchors} as a function of other variables:
    \code{ cbind(Cs, Ce) ~ x1 + x2}}

  \item{data}{a data frame containing two columns Cs, Ce and the covariates
    identified in the formula.}

  \item{weights}{optional case weights in fitting.  Default to 1.}

  \item{start}{initial values for the parameters.  This is in the format
    'c(coefficients, zeta)'}

  \item{...}{additional arguments to be passed to
    \code{\link{optim}[stats]}, most often a 'control' argument.}

  \item{subset}{expression saying which subset of the rows of the data should
    be used in the fit.  All observations are included by default.}

  \item{na.action}{a function to filter missing data.}

  \item{contrasts}{a list of contrasts to be used for some or all of the
    factors appearing as variables in the model formula.}

  \item{Hess}{logical for whether the Hessian (the observed information
          matrix) should be returned. }
  \item{model}{logical for whether the model matrix should be returned. }
  \item{method}{default is probit; alternatives are logistic or
  complementary log-log or cauchit
          (corresponding to a Cauchy latent variable and only available
          in R >= 2.1.0). }
  \item{debug}{additional printing if > 0}
}

\value{An object of classes \code{c("cpolr", "polr")}.  This has
  components

    \item{coefficients}{the coefficients of the linear predictor, which has no
          intercept.}
    \item{zeta}{the intercepts for the class boundaries.}
    \item{deviance}{the residual deviance.}
    \item{fitted.values}{a matrix, with a column for each level of the response.}
    \item{lev}{the names of the response levels.}
    \item{terms}{the 'terms' structure describing the model.}
    \item{df.residual}{the number of residual degrees of freedoms,
      calculated using the weights.}
    \item{edf}{the (effective) number of degrees of freedom used by the
      model.}
    \item{n, nobs}{the (effective) number of observations, calculated using the
          weights. ('nobs' is for use by 'stepAIC').}
    \item{call}{the matched call.}
    \item{convergence}{the convergence code returned by \code{\link{optim}}.}
    \item{niter}{the number of function and gradient evaluations used by
      \code{\link{optim}}.}
    \item{Hessian}{Hessian matrix from \code{\link{optim}}.}

}
\details{

  For cpolr, \code{cpolr.method} default is probit; for additional
	options, see method option in \code{polr}

}

\author{Based on \code{\link[MASS]{polr}} function written by Brian Ripley,
  modifications by Jonathan Wand }

\references{
   Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied
    Statistics with S.}  4th edition.  Springer.

Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., & Wand, J. (2007). Comparing Incomparable Survey Responses: Evaluating and Selecting Anchoring Vignettes. \emph{Political Analysis}, 15(1), 46-66. \doi{10.1093/pan/mpl011}

}
\seealso{\code{\link{anchors}}, \code{\link[MASS]{polr}}}
\examples{

data(freedom)

## an example of directly using cpolr:
ra <- anchors(self ~ vign1 + vign3 + vign6, data = freedom, method ="C")
freedom2 <- insert(freedom, ra )
out <- cpolr(cbind(Cs, Ce) ~ as.factor(country) + sex + educ,
            data = freedom2)
summary(out)


## simplified in the context of anchors:
fo <- list(self= self ~ 1,
           vign = cbind(vign1,vign3,vign6) ~ 1,
           cpolr= ~ as.factor(country) + sex + educ)
ra2 <- anchors(self ~ vign1 + vign3 + vign6, data = freedom, method ="C")
summary(ra, ties="cpolr")

## AVERAGE fitted values
## conditional on observed
fitted(ra2, ties="cpolr", unconditional=FALSE,average=TRUE)
## unconditional prediction
fitted(ra2, ties="cpolr", unconditional=TRUE,average=TRUE)

## fitted probability for each observation
## conditional on observed
fitted(ra2, ties="cpolr", unconditional=TRUE, average=FALSE)
## unconditional prediction
fitted(ra2, ties="cpolr", unconditional=TRUE, average=FALSE)

}

\keyword{models}
