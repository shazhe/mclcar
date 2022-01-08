\name{mcl.HCAR}
\alias{mcl.HCAR}

\title{The Monte Carlo likelihood function of the HCAR model.
}
\description{Evaluate the Monte Carlo likelihood of a hierarchical conditional-autoregressive model at a given parameter value and importance sampler value.}

\usage{
mcl.HCAR(pars, mcdata, data, Evar = FALSE)
}

\arguments{
\item{pars}{parameter values at which the Monte Carlo likelihood function to be evaluated}
  \item{mcdata}{The Monte Carlo sample object returned by the function \code{\link{sim.HCAR}}}
  \item{data}{A list or an environment containing the same variables same as the input for \code{\link{sim.HCAR},}}
  \item{Evar}{If true, the Monte Carlo variance is evaluated.}
}

\value{
\item{mc.lr}{the Monte Carlo likelihood}
\item{v.lr}{if Evar = TRUE, the Monte Carlo variance of the likelihood is also returned}
}

\author{Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{sim.HCAR}}, \code{\link{OptimMCL.HCAR}}, \code{\link{summary.OptimMCL.HCAR}}
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{CAR}
\concept{HCAR}
\concept{Monte Carlo likelihood}
