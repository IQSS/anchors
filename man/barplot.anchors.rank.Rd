\name{barplot.anchors.rank}
\alias{barplot.anchors.rank}
\title{Barlot of distribution of non-parametric ranks}
\description{Barplot of distribution of non-parametric ranks}
\usage{
  \method{barplot}{anchors.rank}(height, ...,
                                 ties = c("uniform","minentropy","omit","cpolr"))
}
\arguments{
   \item{height}{anchors.rank object}
   \item{...}{
     may include additional anchors.rank objects that will plotted
     beside x, but these MUST be listed in consecutive order immediately
     following x.

     in addition, may list extra arguments for \code{barplot} function}

   \item{ties}{shows distribution of ranks by method for
     allocating ties/interval valued cases to a scalar value

     "uniform": allocates in equal proportions among scalar values

     "omit": drops cases with interval values

     "cpolr": uses censored ordered probit model to allocate ties

     "minentropy": allocates ties to produce minimum entropy

   }
}

\author{Jonathan Wand}

\references{
Wand, J., King, G., & Lau, O. (2011). anchors: Software for Anchoring Vignette Data. \emph{Journal of Statistical Software}, 42(3), 1-25. \doi{10.18637/jss.v042.i03}
}
\seealso{\code{\link{anchors}},\code{\link{plot.anchors.rank}}}
\examples{
## see demo(anchors.freedom3)  and demo(anchors.plot) for examples
}
\keyword{hplot}
