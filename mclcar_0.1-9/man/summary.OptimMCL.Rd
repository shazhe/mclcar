\name{summary.OptimMCL}
\alias{summary.OptimMCL}

\title{Summary the output from the iterative procedure of maximising the
  Monte Carlo likelihood.
}
\description{This function summarizes the output of the output from the
  function \code{\link{OptimMCL}}.
}
\usage{
\method{summary}{OptimMCL}(object, family, trace.all = TRUE, mc.covar =
  TRUE, \dots)
}

\arguments{
  \item{object}{ an OptimMCL object returned by \code{\link{OptimMCL}}.}
  \item{family}{ a character takes value in "guass", "binom" and "poisson".}
  \item{trace.all}{ an logic value tells whether the input object given
  by \code{\link{OptimMCL}} contains results from all iterations of not}
  \item{mc.covar}{ if TRUE, the estimated covariance matrix of the
	MC-MLE is returned}
  \item{\dots}{arguments passed to or from other methods.}
}

\value{
  A list containing the following objects:
  \describe{
	\item{MC.mle,}{ the final MC-MLE}
	\item{N.iter,}{ the total number of iterations}
	\item{total.time,}{ the total time elapsed}
	\item{convergence,}{ if TRUE the procedure converges}
	\item{hessian,}{ the Hessian at the MC-MLE if given; the default is NULL}
	\item{mc.covar}{ the estimated covariance matrix of the MC-MLE if
  given; the default is NULL}
	\item{mc.samples}{the Monte Carlo samples size used in the initial
  stage and after the first convergence.}
  }
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{OptimMCL}}
}
\examples{
## See examples for OptimMCL
}

\keyword{CAR}
\keyword{Monte Carlo likelihood}
\keyword{Spatial}
