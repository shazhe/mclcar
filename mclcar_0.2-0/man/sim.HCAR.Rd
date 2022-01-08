\name{sim.HCAR}
\alias{sim.HCAR}

\title{Simulate samples from a HCAR model.
}
\description{Simulate samples from a hierarchical conditional-autoregressive model.}

\usage{
sim.HCAR(psi, data, n.samples)
}

\arguments{
  \item{psi}{A vector of the length of all parameters used for the importance sampler. When the spatial weight matrix W for individual unit is supplied in the data, the first four entries of psi are  
  \describe{
  \item{psi[1] = rho,}{the spatial coefficient parameter for individual units,}
  \item{psi[2] = lambda,}{the spatial coefficient parameter for district level units,}
  \item{psi[3] = sigma.e}{the marginal variance parameter for individual units,}
  \item{psi[4] = sigma.u}{the marginal variance parameter for district level units.}
       }
       Then the rest of entries are values for the linear coefficient beta.
       Otherwise, the function assumes that there is no spatial correlation among the individual unit and the first three entries are used to specify lambda, sigma.e, sigma.u and the rest for beta.
       }
  \item{data}{A list or an environment contains the following variables
  \describe{
  \item{y,}{the observed linear response,}
  \item{X,}{the design matrix for the fixed effects,}
  \item{W,}{the spatial weight matrix for individual units,}
  \item{M,}{the spatial weight matrix for district level units,}
  \item{Z,}{the design matrix for the random effects,}
  \item{In}{an n by n identity matrix, where n is the number of individual spatial units,}
  \item{Ik}{a k by a indentity matrix, where k is the number of distrct level unites.}
    }}
  \item{n.samples}{the number of Monte Carlo samples}
}

\value{
\item{psi}{the importance sampler value used for generating the samples,}
\item{n.samples}{the number of Monte Carlo samples,}
\item{u.y}{the Monte Carlo samples of the randome effects,}
\item{lpsi}{the unnormalsied likelihood components evaluated at psi,}
\item{const.psi}{the constant in the likelihood function evaluated at psi}
}

\author{
Zhe Sha \email{zhesha1006@gmail.com}
}

\seealso{
\code{\link{mcl.HCAR}}, \code{\link{OptimMCL.HCAR}}, \code{\link{summary.OptimMCL.HCAR}}
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{CAR}
\concept{HCAR}
\concept{Monte Carlo likelihood}
