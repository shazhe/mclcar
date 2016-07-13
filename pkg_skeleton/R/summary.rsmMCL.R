#### summary function for rsmMCL
summary.rsmMCL.lm <- function(object, trace.all = TRUE, mc.covar=TRUE){
    niter <- object$N.iter
    if(trace.all){
        rsm.mle <- object$Psi[[niter]]
        mc.data <- object$Sim.data[[niter]]
    }else{
        rsm.mle <- object$Psi[[2]]
        mc.data <- object$Sim.data[[1]]
    }

    total.time <- object$total.time
    convergence <- object$convergence
    mc.samples <- object$mcsamples
    data <- object$data

    if(mc.covar){
        B <- Hessian.dCAR(rsm.mle, object$data) # Hessian at the mc.mle
        invB <- solve(B)
        A <- MCMLE.var(rsm.mle, data = data, simdata = mc.data)
        mc.var <- invB %*% A %*% invB/nrow(data$data.vec)
        ans <- list(MC.mle = rsm.mle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = B, mc.covar = mc.var,
                    mc.samples = mc.samples)
    }else{
        ans <- list(MC.mle = rsm.mle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = NULL, mc.covar = NULL,
                    mc.samples = mc.samples)
    }
    return(ans)
}

summary.rsmMCL.glm <- function(object, family, trace.all = TRUE, mc.covar=TRUE){
    niter <- object$N.iter
    if(trace.all){
        rsm.mle <- object$Psi[[niter]]
        mc.data <- object$Sim.data[[niter]]
        mc.var <- object$MC.Vars[[niter]]

    }else{
        rsm.mle <- object$Psi[[2]]
        mc.data <- object$Sim.data[[1]]
        mc.var <- object$MC.Vars[[1]]
    }

    total.time <- object$total.time
    convergence <- object$convergence
    mc.samples <- object$mcsamples
    ##data <- object$data


    if(mc.covar){
        Hess <- mcl.Hessian.glm(rsm.mle, mcdata = mc.data, family = family)
        opt.res <- list(estimate = rsm.mle, hessian = Hess)
        mc.var <- vmle.glm(opt.res, mcdata = mc.data, family = family)
        ans <- list(MC.mle = rsm.mle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = Hess, mc.covar = mc.var,
                    mc.samples = mc.samples)
    }else{
        ans <- list(MC.mle = rsm.mle, N.iter = niter, total.time = total.time,
                    convergence = convergence, hessian = NULL, mc.covar = NULL,
                    mc.samples = mc.samples)
    }
    return(ans)
}

summary.rsmMCL <- function(object, family, trace.all = TRUE, mc.covar = TRUE, ...){
    if(family == "gauss"){
        summary.rsmMCL.lm(object, trace.all, mc.covar)
    }else{
        summary.rsmMCL.glm(object, family, trace.all, mc.covar)
    }
}
