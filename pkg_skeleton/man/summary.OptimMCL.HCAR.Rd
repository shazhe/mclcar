\name{summary.OptimMCL.HCAR}
\alias{summary.OptimMCL.HCAR}
\alias{ranef.HCAR}

\title{Summary the output from the iterative procedure of maximising the
  Monte Carlo likelihood.
}
\description{The summary function summarizes the output of the output from the
  function \code{\link{OptimMCL.HCAR}} and the \code{ranef.HCAR} calculate the empirical Bayesian estimates of the random effects given the Monte Carlo maximum likelihood estimates.
}
\usage{
\method{summary}{OptimMCL.HCAR}(object, trace.all = TRUE, mc.covar =
  TRUE, \dots)
  
ranef.HCAR <-function(pars, data)
}

\arguments{
  \item{object}{ an OptimMCL object returned by \code{\link{OptimMCL.HCAR}}.}
  \item{trace.all}{ an logic value tells whether the input object given
  by \code{\link{OptimMCL.HCAR}} contains results from all iterations of not}
  \item{mc.covar}{ if TRUE, the estimated covariance matrix of the
	MC-MLE is returned}
  \item{\dots}{arguments passed to or from other methods.}
  
  \item{pars}{the paramter values for calculating the empirical Bayesian estimates of the randome effects.}
  \item{pars}{a list or enivironment of data same as described in \code{\link{sim.HCAR}}}
}

\value{
  The summary function returns a list containing the following objects:
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
  
  The \code{ranef.HCAR} function return a dataframe object containing the estimated random effects and their corresponding standard deviations.
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{mcl.HCAR}}, \code{\link{sim.HCAR}}, \code{\link{OptimMCL.HCAR}}
}
\examples{
## See examples for OptimMCL
}

\keyword{HCAR}
\keyword{Monte Carlo likelihood}
\keyword{Spatial}
