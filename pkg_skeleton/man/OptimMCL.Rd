\name{OptimMCL}
\alias{OptimMCL}
\title{Iterative procedure for maximising the Monte Carlo likelihood}
\description{
The function uses an iterative procedure of directly maximising the
Monte Carlo likelihood and the updating step size is limited by defining
an experimental region using the estimated Monte Carlo variance.
}
\usage{
OptimMCL(data, psi0, family, control = list(), mc.control = list())
}

\arguments{
  \item{data}{ the data set to be estimated; for direct CAR model it is
  a list objects same used in \code{\link{loglik.dCAR}} and for glm with CAR
  latent variables, it is the same as the output of \code{\link{CAR.simGLM}}.
}
  \item{psi0}{ the initial value for the importance sampler parameter}
  \item{family}{a character take values in "gauss" (for direct CAR
  models), "binom" and "poisson" (glm Binomial/Poisson with CAR
  latent variable).}
 \item{control}{ a list that controls the iterative procedures, see more
  in Details.}
  \item{mc.control}{ a list that controls the MCMC algorithm, see more
  in Details of \code{\link{postZ}}}
}

\details{
  The iterative procedure starts by using samples from the importance
  sampler with the parameter values given by \code{psi0} and directly
  optimise the Monte Carlo likelihood. Then the importance sampler
  parameter value is updated by using either the MLE found by current
  step if it is within the experimental region (to be defined in the
  following), or a weighted average between the current psi value and
  the current MLE value, with weight to be specified by \code{psi.ab} in
  the list of \code{control} argument.

  The experimental region is defined to be the region where the
  estimated variance of the Monte Carlo likelihood is smaller than
  \eqn{a/\sqrt(n.s)} where \eqn{n.s} is the number of Monte Carlo
  samples and \eqn{a} and is some fixed constant that can be specified
  through \code{psi.ab} in \code{control}.

  The procedure iterate until convergence which is defined to be that
  the Monte Carlo likelihood is smaller than min(1, 2*sqrt(mcl.var))
  where mcl.var is the variance of the estimated variance of the Monte
  Carlo likelihood evaluated at the MC-MLE. At this point, the procedure
  use a increased number of Monte Carlo samples and do the same
  iteration until converges again. The MC-MLE found on this second
  convergence is given as the final result. The Monte Carlo sample size
  is increased by a multiple of user specified number through the
  \code{control} argument.

  The \code{control} argument is a list containing the following
  elements:
  \describe{
	\item{n.iter,}{ maximum number of iterations; default value is 10
  for direct CAR models and 20 for glm with latent CAR variables}
	\item{mc.samples,}{ (for direct CAR model only) an array of two with
  the first element being the number of Monte Carlo samples and the
  second being the increasing multiple of the sample size after the first
  convergence; default is c(500, 2)}
\item{s.increase,}{ (for glm with latent CAR variables only) an numeric
 value for the increasing multiple of the sample size after the first
 convergence; default value is 2.}
\item{rho.range,}{ the valid range for the parameter rho in the direct CAR
  model; default value is c(-0.249, 0.249)}
\item{psi.range,}{ the valid range for the parameters rho and sigma in
  the glm with latent CAR variables; default is c(-0.249, 0.249, 0.1, 10)}
\item{psi.ab,}{ an array with two entries controlling the experimental
  region and updating rules; the first entry controls the experimental
  region as \eqn{a} defined above and the second entry controls the
  updating rule as defined above; default is c(1, 0.5)}
\item{mc.var,}{ if true, the estimated covariance matrix for the MC-MLE
  is returned in every iteration; default is FALSE}
\item{trace.all,}{ if true, each iteration is stored and returned in the
  final output, default is TRUE}
\item{verbose,}{ if true, an summary message is printed after each
  iteration, default is TRUE}
}
}

\value{When trace.all is TRUE, the function returns a list containing the following objects:
  \describe{
	\item{Psi,}{ a list of importance sampler parameter values from all iterations;}
	\item{MC.datas,}{ a list object generated by
  \code{\link{mcl.prep.dCAR}} or \code{\link{mcl.prep.glm}} from all iterations;}
	\item{MC.MLEs,}{ a list of the MC-MLE found by all iterations;}
	\item{MC.Hess,}{ a list of the Hessian matrix at the MC-MLE from all
  iterations;}
	\item{MC.Vars,}{ a list of the estimated covariance matrix of the
  MC-MLE from all iterations;}
	\item{data,}{ the data object supplied to the function;}
	\item{N.iter,}{ the total number of iterations;}
	\item{total.time,}{ the total time elapsed;}
	\item{convergence,}{ a logical value indicating whether the
  procedure converged or not;}
	\item{mcsamples,}{ an array of two entries for the initial Monte
  Carlo sample size and the increased Monte Carlo sample size after the
  first convergence.}
	}

	When trace.all is FALSE, the function returns the same list with
	each object in the list containing the result of the final iteration only.
  }

\references{Sha, Z. 2016 \emph{Estimating conditional auto-regression models}, DPhil
Thesis, Oxford.
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{summary.OptimMCL}}, \code{\link{rsmMCL}}
}
\examples{

## Take long time to run
## Simulate some data to work with
set.seed(33)
n.torus <- 10
nb <- 30
rho <- 0.2
sigma <- 1.5
beta <- c(1, 1)
pars.true <- c(rho, sigma, beta)
X0 <- cbind(rep(1, n.torus^2), sample(log(1:n.torus^2)/5))
mydata2 <- CAR.simGLM(method = "binom", n = c(n.torus, n.torus), pars = pars.true,
                      Xs = as.matrix(X0), n.trial = nb)

## use a glm to find initial values for the importance sampler
data.glm <- data.frame(y=mydata2$y, mydata2$covX[,-1])
fit.glm <- glm(cbind(y, nb-mydata2$y) ~ .,data = data.glm, family=binomial)
library(spatialreg)
library(spdep)
logitp <- log((mydata2$y+0.5)/(mydata2$n.trial - mydata2$y + 0.5))
data.splm <- data.frame(y=logitp, mydata2$covX[,-1])
listW <- mat2listw(mydata2$W)
fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
pars1 <- c(fit.splm$lambda, fit.splm$fit$s2, coef(fit.glm))

## Use the iterative procedure to find the MC-MLE
## !!!NOTE: the example below is only an illustration of usage
## users should increase the number of iterations and MCMC samples
## to get convergence results.
iter.mcmle <- OptimMCL(data = mydata2, psi0 = pars1, family = "binom",
                       control = list(n.iter = 1, mc.var = TRUE),
                       mc.control = list(N.Zy = 1e3, Scale = 1.65/(n.torus^(2/6)), thin = 5,
                   burns = 5e2, method = "mala", scale.fixed = TRUE))
summary(iter.mcmle, family = "binom", mc.covar=TRUE)
}

\concept{Monte Carlo likelihood}
\concept{CAR}
\concept{Spatial}