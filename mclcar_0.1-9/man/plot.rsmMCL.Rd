\name{plot.rsmMCL}
\alias{plot.rsmMCL}
\title{
Plot the fitted response surfaces
}
\description{
This function plot the fitted response surfaces by using the result from \code{\link{rsmMCL}}.
}
\usage{
\method{plot}{rsmMCL}(x, family, trace.all = TRUE, plain=TRUE, \dots)
}
\arguments{
  \item{x}{ an rsmMCL object retruned by \code{\link{rsmMCL}}.
}
  \item{family}{
a character takes value in "guass", "binom" and "poisson".
}
  \item{trace.all}{
an logic value tells whether the input object given
by \code{\link{rsmMCL}} contains results from all iterations or not.
}
\item{plain}{
  an logic value indicates whether to show plots in separate devices
  with a smart display. Use the defalut value TRUE if the function
  is used in generating documents by using \pkg{knitr}.
}
\item{\dots}{For extra arguments. The Lattice parameter \code{layout}
  can be used to control the layout of the graphs.}

}


\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{rsmMCL}}, \code{\link{summary.rsmMCL}}
}
\examples{
## See examples for rsmMCL
}
\keyword{CAR}
\keyword{Monte Carlo likelihood}
\keyword{Spatial}
\keyword{Response surface design}