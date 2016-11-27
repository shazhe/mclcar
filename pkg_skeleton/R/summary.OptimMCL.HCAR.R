#### summary function for OptimMCL
summary.OptimMCL.HCAR <- function(object, trace.all = TRUE, mc.covar=TRUE){
    niter <- object$N.iter
    if(trace.all){
        mcmle<- object$MC.MLEs[[niter]]
        mle.hessian <- object$MC.Hess[[niter]]
        mc.var <- object$MC.Vars[[niter]]
    }else{
        mcmle <- object$MC.MLEs
        mle.hessian <- object$MC.Hess
        mc.var <- object$mc.Vars
    }

    total.time <- object$total.time
    convergence <- object$convergence
    mc.samples <- object$mcsamples
    ##data <- object$data


    if(mc.covar){
        ans <- list(MC.mle = mcmle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = mle.hessian, mc.covar = mc.var,
                    mc.samples = mc.samples)
    }else{
        ans <- list(MC.mle = mcmle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = NULL, mc.covar = NULL,
                    mc.samples = mc.samples)
    }
    return(ans)
}


## Random effect estimate
ranef.HCAR <-function(pars, data){
  if (is.environment(data)) {
    y <- get("y", envir=data)
    X <- get("X", envir=data)
    W <- get("W", envir=data)
    M <- get("M", envir=data)
    Z <- get("Z", envir=data)
    n <- length(y)
    K <- nrow(M)
    In <- get("In", envir=data) #diag(1, n)
    Ik <- get("Ik", envir=data) #diag(1, K)
   } else {
    y <- data$y
    X <- data$X
    W <- data$W
    M <- data$M
    Z <- data$Z
    n <- length(y)
    K <- nrow(M)
    In <- data$In #diag(1, n)
    Ik <- data$Ik #diag(1, K)
  } 
 

    if(is.null(W)){
    lambda <- pars[1]
    sigma.e <- pars[2]
    sigma.u <- pars[3]
    beta <- pars[-c(1:3)]
    A <- In
    B <- Ik - lambda*M
}else{
    rho <- pars[1]
    lambda <- pars[2]
    sigma.e <- pars[3]
    sigma.u <- pars[4]
    beta <- pars[-c(1:4)]
    A <- In - rho*W
    B <- Ik - lambda*M
}

    Q.u <- B/sigma.u
    Q.e <-  A/sigma.e

    res <- y - X %*% beta
    Cov.b <- solve(t(Z) %*% Q.e %*% Z + Q.u)

    b <- Cov.b %*% t(Z) %*% Q.e %*% res
    sd.b <- sqrt(diag(Cov.b))
    return(data.frame(ranef = b, ranef.sd = sd.b))
    }
