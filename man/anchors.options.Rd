\name{anchors.options}
\alias{anchors.options}
\title{Set or query anchors() parameters}

\description{ 'anchors.options' can be used to set or query parameters passed to
  anchors functions.  Parameters can be set by specifying them as
  arguments to 'anchors.options' in 'tag = value' form, and will return an
  updated list of parameters.
}

\usage{
  anchors.options(...)
}

\arguments{
  \item{...}{arguments in 'tag = value' form.
    The tags must come from the anchors parameters described
    below.}

}
\section{anchors.options()}{

  Options specific to each method are described in their respective help
  pages.

  Options generic to all functions are:

    \describe{
      \item{verbose     }{ FALSE, if TRUE more detailed run-type printing}
      \item{silence     }{ FALSE, if TRUE turns off even some standard printing}
      \item{debug       }{ 0, if >0 show gory details for
	debugging. Caution: high numbers can produce very high volume output.}
    }

}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}

\seealso{
  \code{\link{anchors}},
  \code{\link{chopit}},
  \code{\link{anchors.data}},
  \code{\link{anchors.order}}
}

\examples{

data(mexchn)
fo <- list(self = xsayself ~ 1,
           vign = xsay3    ~ 1,
           tau  =           ~1 )
out1  <- chopit( fo, mexchn, options= anchors.options( vign.var = "hetero", normalize="hilo" ))
summary(out1)


}

\keyword{models}
