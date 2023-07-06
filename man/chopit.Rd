\name{chopit}
\alias{chopit}
\title{Compound Hierarchical Ordered Probit (CHOPIT)}
\description{
Estimate Compound Hierarchical Ordered Probit (CHOPIT) model
}
\usage{
  chopit(formula, data, subset,
         options=anchors.options(),
         na.action = na.omit)
}
\arguments{
  \item{formula}{A list of named formulas giving a symbolic description
  of the model to be fit. See 'Details' below.}
  \item{data}{An anchors.data class object.  See 'Details' below.}
  \item{subset}{Logical expression indicating elements or rows to keep:
          missing values are taken as false; equivalent to function \code{subset}.}
  \item{options}{An anchors.options class object.  See 'Details'
  below.}
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

      Only options relevant to \code{chopit} are discussed here.

    The named list() of formulas may include

    self: self-response variable and linear predictors.

    vign: vignette responses, bound together by cbind (as in glm
    binomial syntax), with ONLY '~ 1'  as RHS of equation (formulas must
    have RHS).  Must be present for all methods.

    tau : linear predictors of cutpoints, NO LHS variable.

    tau1 : Optional.  linear predictors of only the mean shift of all cutpoints
    NO LHS variable. Default: if omitted, set equal to
    formula specified in 'tau='.

    Example:

    \preformatted{
    fo <- list(self = xsayself ~ age,
              vign = cbind(xsay3,xsay1)    ~ 1,
              tau  = ~ age + china,
              tau1 = ~ age + china + male  )
	      }

    One can also use a single formula as a short hand: specifying

    fo <- self ~ vign1 + vign2

    is equivalent to

    \preformatted{
    fo <- list(self = self ~ 1
               vign = cbind(vign1,vign2) ~ 1,
               tau  = ~ 1,
               tau1 = ~ 1,
               cpolr= ~ 1)
	     }


  All of the response variables must be in the form of consequetive non-negative
  numeric integers, i.e., 1, 2, ... K.

    ** \code{anchors} currently does not support factor responses.**

  The \code{chopit} model does *not* by default require that every
  respondent answer all vignette or the self questions to be included in
  the model; only answers that are non-missing are used in the
  likelihood function.

}


\section{options and anchors.options()}{

  Invoking the chopit model without options,

  \code{ chopit(formula, data)}

  is equivalent to,

  \preformatted{chopit(formula, data,
                  options = anchors.options( normalize ="self",
                                             vign.var  = "hetero",
                                             linear    = TRUE,
                                             analytical= TRUE,
                                             optimizer= "optim"))
					   }



 To change an option setting, pass one or more of the following tagged
 values to the options argument via the \code{anchors.options()}
 function.  Alternative vluaes for options are passed using the
 tag=new.value syntax (same with the \code{par()} function.

 To see all default values, type \code{anchors.options()}
 without arguments.


 \describe{
  \item{normalize}{Normalization constraints
    imposed.

    "self" = location is set by removing any intercept from x0
    vector and constraining variance of first self question to be 1.

    "hilo" = sets location and scale by setting first theta to 0 and last
    theta to 1.
  }

  \item{vign.var}{
    "homo" constrains all vignettes to have a single, common variance.

    "hetero" allows variances to differ for each vignette}


  \item{linear}{
    TRUE: (default) use additive linear formulation of taus if true.

    FALSE: uses additive exponentiated values to force each cutpoint to be
    positive.

    NOTE: Analytical gradients only availabe with linear formulation.}

  \item{analytical}{

    TRUE; (default) use analytical gradients--much faster than
    numerical.

    FALSE: use numerical gradients

    Currently analytical gradients are only available for model with
    normalize='self' and linear=TRUE options.
  }

  \item{optimizer   }{

    "optim": (default) estimate using \code{optim}

    "genoud": estimate chopit using \code{\link[rgenoud]{genoud}}
    optimizer

  }
  }

  With optimizer = "optim", there are additional options for specifying
  how \code{\link{optim}} is invoked.
  \describe{
    \item{optim.method}{This sets the 'method' of \code{\link{optim}}
      and can be one of

      optim.method = c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN")

      "BFGS" is the default.

      The Details section of  \code{\link{optim}} gives the precise
      definition of each method.}

    \item{maxit       }{This option controls the maximum number of
      interations that optim will perform before stopping, even if a
      solution has not been found.  Default: 500.}
  }

  With optimizer = "genoud", there are additional options for specifying
  how \code{\link[rgenoud]{genoud}} is invoked.
  \describe{
    \item{wait.generations}{If there is no improvement in the objective function
          in this number of generations, 'genoud' will think that it
          has found the optimum.  Default: 1}
    \item{pop.size    }{ Population Size.  This is the number of individuals 'genoud'
      uses to solve the optimization problem. Default: 500}

    \item{MemoryMatrix}{ This variable controls if 'genoud' sets up a memory
          matrix.  Such a matrix ensures that 'genoud' will request the
          fitness evaluation of a given set of parameters only once.
          The variable may be 'TRUE' or 'FALSE'.  If it is 'FALSE',
          'genoud' will be aggressive in conserving memory. Default:
	  TRUE. }
    \item{max.generations }{ Maximum Generations.  This is the maximum number of
          generations that 'genoud' will run when attempting to
          optimize a function.  Default: 100. }
    \item{domain      }{'Genoud' will create a Domains
          matrix by setting the lower bound for all of the parameters
          equal to -1 * 'default.domains' and the upper bound equal to
          'default.domains'.  Default:  5,}

      }


}


\value{

  \item{data}{list of class anchors.data}
  \item{parm}{list which contains MLEs }
  \item{count}{list containing summary of data dimensions}
  \item{options}{list of options, possibly modified from original call
    to anchors() or internal consistency with use of chopit}
  \item{optim}{complete object returned by optim() or genoud() call}
  \item{hess}{Hessian}
  \item{LL.vign}{likelihood values for the vignette component of model}
  \item{LL.self}{likelihood values fo rthe self component of model}
  \item{gr}{gradients of fitted model}
  \item{time}{timing information for estimation}
}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}

King, G., Murray, C. J. L., Salomon, J., & Tandon, A.. (2004). Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research. \emph{American Political Science Review}, 98(1), 191-207. \doi{10.1017/S000305540400108X}
}


\seealso{\code{\link{anchors}},
  \code{\link{anchors.data}},
  \code{\link{anchors.options}},
  \code{\link{optim}},
  \code{\link[rgenoud]{genoud}}}
\examples{

##
data(mexchn)
fo <- list(self = xsayself ~ 1,
           vign = xsay3    ~ 1,
           tau  =           ~1 )
out1  <- chopit( fo, mexchn)
summary(out1)


}
\keyword{models}
