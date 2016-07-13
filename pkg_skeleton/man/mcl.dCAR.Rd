\name{mcl.dCAR}
\alias{mcl.dCAR}
\alias{mcl.prep.dCAR}
\alias{mcl.profile.dCAR}
\alias{vmle.dCAR}
\alias{Avar.lik.dCAR}

\title{
Monte Carlo likelihood calculation for direct CAR models.
}
\description{This help page contains functions for calculating the Monte
  Carlo likelihoods and variances for direct CAR models.

  \code{mcl.prep.dCAR} generates the Monte Carlo samples from the
  importance sampling distribution to be used in the Monte Carlo
  likelihood functions.

  \code{mcl.dCAR} calculates the Monte Carlo log-likelihood ratio for a given parameter
  value to the importance sampler value.

  \code{mcl.profile.dCAR} calculates the Monte Carlo profile
  log-likelihood ratio of rho.

  \code{vmle.dCAR} calculates the variance at the Monte Carlo MLE.

  \code{Avar.lik.dCAR} calculates the asymptotic variance of the Monte
  Carlo likelihood for a direct CAR model on a given size of torus and
  given size of Monte Carlo samples.
}
\usage{
mcl.prep.dCAR(psi, n.samples, data)

mcl.dCAR(pars, data, simdata, Evar = FALSE)
mcl.profile.dCAR(rho, data, simdata, Evar = FALSE)

vmle.dCAR(MLE, data, simdata)
Avar.lik.dCAR(pars, psi, data, n.samples, Log = TRUE)
}

\arguments{
  \item{psi}{the value of the importance sampler parameters}
  \item{n.samples}{the number of Monte Carlo samples}
   \item{pars}{the parameter value for the Monte Carlo likelihood function to be
	 evaluated at}
   \item{rho}{the value of rho for the Monte Carlo profile likelihood
  function to be evaluated at}
  \item{data}{a \code{data} object same as described in
	\code{\link{loglik.dCAR}}}
  \item{simdata}{a list object in the same format as the output of
	\code{mcl.prep.dCAR}}
  \item{Evar}{if TRUE return the estimated variance of the
	Monte Carlo likelihood}
  \item{MLE}{a list object with
	\describe{
	  \item{\code{par}}{the estimated Monte Carlo MLE}
	  \item{\code{hessian}}{the hessian at the Monte Carlo MLE}
	}
	}
  \item{Log}{if log = TRUE, then the variance is of the
	log-likelihood ratio; otherwise it is of the likelihood ratio}
}

\details{The asymptotic and estimated variance is derived for a direct CAR model on a \eqn{n
	\times n} torus with `rook` style binary neighbourhood
  matrix. Details of the derivation is given in Sha (2016) ch. 3.
}
\value{\code{mcl.prep.dCAR} returns a list object containing the
  following following elements:
  \itemize{
	\item{\code{simY}}{a matrix of the Monte Carlo samples}
    \item{\code{t10} and \code{t20}}{statistics pre-calculate for the importance
	  distribution}
  }

  \code{mcl.dCAR} and \code{mcl.proifle.dCAR} return a numeric value of the Monte Carlo likelihood
  when \code{Evar = FALSE}; if \code{TRUE} it returns an array of the
  Monte Carlo likelihood and the corresponding estimated variance.

  \code{vmle.dCAR} returns the estimated covariance matrix of the Monte
  Carlo MLE.

  \code{Avar.lik.dCAR} returns a numeric value of the asymptotic
  variance.
}
\references{
Sha, Z. 2016 \emph{Estimating conditional auto-regression models}, DPhil
Thesis, Oxford.
}
\author{
Zhe Sha \email{sha@stats.ox.ac.uk}
}


\seealso{
\code{\link{mcl.glm}}, \code{\link{OptimMCL}}, \code{\link{rsmMCL}}
}
\examples{
## Simulate some data to work with
set.seed(30)
n.torus <- 10
rho <- 0.15
sigma <- 1.5
beta <- c(1, 2)
XX <- cbind(rep(1, n.torus^2), log(1:n.torus^2))
mydata1 <- CAR.simTorus(n.torus, n.torus, rho, 1/sigma)
y <- XX \%*\% beta + mydata1$X
mydata1$data.vec <- data.frame(y=y, XX[,-1])

## Choose the importance sampler value to the the MPLE
psi1 <- mple.dCAR(mydata1)

## Prepare the Monte Carlo samples
mcdata1 <- mcl.prep.dCAR(psi = psi1, n.samples = 100, data = mydata1)

## Calculate the Monte Carlo likelihoods
mcl.dCAR(c(rho, sigma, beta), data = mydata1, simdata = mcdata1, Evar =
TRUE)
mcl.profile.dCAR(rho, data = mydata1, simdata = mcdata1, Evar = TRUE)

## Do a direct optimization of the Monte Carlo likelihood function
## to find the MLE but it takes very long time to run and may not converge
\dontrun{
opt <- optim(psi1, fn = mcl.dCAR,
             data = mydata1, simdata = mcdata1,
             hessian = TRUE, control = list(fnscale = -0.5))
}

## Assume the we have obtained the opt with the following values
opt<- list(par = c (0.08013547, 1.70294099, 0.01571957, 2.23203089),
           hessian =  matrix(c(  -190.8791352,  -1.5682773,    0.1863733,   -4.844151,
 -1.5682773,  -16.1605186,    0.4911009,    4.403844,
 0.1863733,   0.4911009,  -41.1496774, -156.686631,
-4.8441507,   4.4038442, -156.6866312, -634.316017), 4, 4))

## Calculate the variance of the MC-MLE
vmle.dCAR(opt, data = mydata1, simdata = mcdata1)

## Calculate the asymptotic variance of the likelihood at MC-MLE
Avar.lik.dCAR(pars = opt$par, psi = psi1, data = mydata1, n.samples =
100)

}

\keyword{Monte Carlo variance}
\keyword{Monte Carlo likelihood}
\keyword{CAR}