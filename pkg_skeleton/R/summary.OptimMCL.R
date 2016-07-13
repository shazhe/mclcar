#### summary function for OptimMCL
summary.OptimMCL.lm <- function(object, trace.all = TRUE, mc.covar=TRUE){
    niter <- object$N.iter
    if(trace.all){
        mcmle.rho <- object$MC.MLEs[[niter]]
        mle.hessian <- object$MC.Hess[[niter]]
        mc.datas <- object$MC.datas[[niter]]
    }else{
        mcmle.rho <- object$MC.MLEs
        mle.hessian <- object$MC.Hess
        mc.datas <- object$MC.datas
    }

    total.time <- object$total.time
    convergence <- object$convergence
    mc.samples <- object$mcsamples
    data <- object$data
    mcmle.bs <- sigmabeta(mcmle.rho, object$data)
    mcmle <- c(mcmle.rho, mcmle.bs)

    if(mc.covar){
        B <- Hessian.dCAR(mcmle, object$data) # Hessian at the mc.mle
        invB <- solve(B)
        A <- MCMLE.var(mcmle, data = data, simdata = mc.datas)
        mc.var <- invB %*% A %*% invB/nrow(data$data.vec)
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

#### summary function for OptimMCL
summary.OptimMCL.glm <- function(object, trace.all = TRUE, mc.covar=TRUE){
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

#### wrapper for direct and latent models
summary.OptimMCL <- function(object, family, trace.all = TRUE, mc.covar = TRUE, ...){
    if(family == "gauss"){
        summary.OptimMCL.lm(object, trace.all, mc.covar)
    }else{
        summary.OptimMCL.glm(object,trace.all, mc.covar)
    }
}
