################################################################
####   Direct Optimisation using an iterative procedure     ####
################################################################
#### Optim function using profile log-likelihood of rho for direct car models
OptimMCL.lm <- function(data, psi0, control = list()){
    con <- list(mc.samples = c(500,2), rho.range = c(-0.249, 0.249),
                n.iter = 10, psi.ab = c(1,0.5), mc.var = FALSE,
                trace.all = TRUE, verbose = TRUE)
    ##namc <- names(con)
    ##con[(namc <- names(control))] <- control
    con[names(control)] <- control

    ## Initializing
    n.s <- n.s0 <- con$mc.samples[1]
    n.increase <- con$mc.samples[2]
    n.iter = con$n.iter
    psi.a <- con$psi.ab[1]
    psi.b <- con$psi.ab[2]
    psi <- psi0
    rho.range <-con$rho.range
    n.data <- nrow(data$W)

    Psi <- list()
    MC.datas <- list()
    MC.MLEs <- list()
    MC.Hess <- list()
    MC.SDs <- list()
    flag.converge <- FALSE
    i = 1
    t.total <- 0

    while(i <= n.iter & !flag.converge){
        tt <- system.time({
                           mcdata <- mcl.prep.dCAR(psi, n.s, data)
                           opt <- optim(psi[1], mcl.profile.dCAR, data = data, simdata = mcdata,
                                        method = "L-BFGS-B",
                                        lower = rho.range[1], upper = rho.range[2],
                                        hessian = TRUE, control = list(fnscale = -1))
                           mc.mle <- opt$par
                           mc.Hess <- opt$hessian
                           mcl <- mcl.profile.dCAR(mc.mle, simdata = mcdata, data = data, Evar = TRUE)
                           mcl.var <- mcl[2]
                           mle.flag <- ifelse(is.na(mcl.var), FALSE, mcl.var <= psi.a/n.s)
                           sigmabeta.new <- sigmabeta(mc.mle, data = data)
                           if(mle.flag){
                               psi.new <- c(mc.mle, sigmabeta.new)
                               if(con$mc.var){
                                   B <- Hessian.dCAR(psi.new, data) # Hessian at the mc.mle
                                   invB <- solve(B)
                                   A <- MCMLE.var(psi.new, data = data, simdata = mcdata)
                                   mc.var <- invB %*% A %*% invB
                                   mc.sd <- sqrt(mc.var[1,1]/n.data)
                               }else{
                                   mc.sd <-NA
                               }
                           }else{
                               psi.new <- psi.b*(psi + c(mc.mle, sigmabeta.new))
                               mc.sd <- NA
                           }
                       })
        Psi[[i]] <- psi.new
        MC.datas[[i]] <- mcdata
        MC.MLEs[[i]] <- mc.mle
        MC.Hess[[i]] <- mc.Hess
        MC.SDs[[i]] <- mc.sd

        psi <- psi.new
        flag.converge1 <- (mcl[1] < 1) & (mcl[1] < 2*sqrt(mcl.var))
        flag.converge <- flag.converge1 & (n.s > n.s0)
        t.total <- t.total + tt[3]
        if(con$verbose){
            print(paste0(c("i = ", "n.s = ", "mc.mle = ", "converge = ", "time = "),
                         c(i, n.s, psi.new[1], flag.converge, t.total)))
        }
        n.s <- ifelse(flag.converge1, n.s*n.increase, n.s)
        i <- i+1
    }
    if(con$trace.all){
        ans <- list(Psi=Psi, MC.datas=MC.datas, MC.MLEs=MC.MLEs, MC.Hess=MC.Hess,
                    MC.SDs=MC.SDs, data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s/n.increase))
        attr(ans, "class") <- "OptimMCL"
        return(ans)
    }else{
        ans <- list(Psi=Psi[[i-1]], MC.datas=MC.datas[[i-1]], MC.MLEs=MC.MLEs[[i-1]],
                    MC.Hess=MC.Hess[[i-1]], MC.SDs=MC.SDs[[i-1]],
                    data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "OptimMCL"
        return(ans)
    }
}


#### Optim function using profile log-likelihood of rho for latent car models
OptimMCL.glm <- function(data, psi0, family, control = list(), mc.control = list()){
    N.data <- nrow(data$W)
    ## MCMC default control
    mc.con <- list(method = "mala", N.Zy = 1e4, Scale = 1.65/(N.data^(2/6)), thin = 5,
                   burns = 5e3, scale.fixed = TRUE, c.ad = c(1, 0.7))
    ##namc.mc <- names(mc.con)
    ##mc.con[(namc.mc <- names(mc.control))] <- mc.control
    mc.con[names(mc.control)] <- mc.control
    ## general algorithm default control
    con <- list(n.iter = 20, s.increase = 2, psi.ab = c(1,0.5),
                psi.range = c(-0.249, 0.249, 0.1, 10), mc.var = FALSE,
                trace.all = TRUE, verbose = TRUE, print = 0)
    ##namc <- names(con)
	##con[(namc <- names(control))] <- control
	con[names(control)] <- control

    ## Initializing
    n.s <- n.s0 <- (mc.con$N.Zy - mc.con$burns)/mc.con$thin
    n.increase <- con$s.increase
    n.iter = con$n.iter
    psi.a <- con$psi.ab[1]              # constant a in updating criterion 3.6
    psi.b <- con$psi.ab[2]              # constant b in updating criterion 3.6
    psi <- psi0
    psi.range <- con$psi.range

    ## Constraints on parameters
    A <- matrix(0, nrow = 4, ncol = length(psi0))
    A[,1] <- c(1, -1, 0, 0)             # for rho
    A[,2] <- c(0, 0, 1, -1)             # for sigma
    B <- c(-psi.range[1], psi.range[2], -psi.range[3], psi.range[4])

    Psi <- list()
    MC.datas <- list()
    MC.MLEs <- list()
    MC.Hess <- list()
    MC.Vars <- list()
    flag.converge <- FALSE
    i = 1
    t.total <- 0

    while(i <= n.iter & !flag.converge){
        tt <- system.time({mc.data <- mcl.prep.glm(data = data, family = family,
                                                  psi = psi, mcmc.control = mc.con)
                           opt.res <- maxBFGS(fn = mcl.glm, grad = mcl.grad.glm,
                                              hess = mcl.Hessian.glm,
                                              family = family, mcdata = mc.data,
                                              start=as.numeric(psi), print.level=con$print,
                                              constraints=list(ineqA=A, ineqB=B))
                           mc.mle <- opt.res$estimate
                           mc.Hess <- opt.res$hessian
                           mcl <- mcl.glm(mc.mle, mcdata = mc.data, family = family, Evar = TRUE)
                           mcl.var <- mcl[2]
                           mle.flag <- ifelse(is.na(mcl.var), FALSE, round(mcl.var*n.s) <= psi.a)
                           if(mle.flag){
                               psi.new <- mc.mle
                               if(con$mc.var){
                                   mc.var <- vmle.glm(opt.res, mcdata = mc.data, family = family)
                               } else{
                                   mc.var <-NA
                               }
                           } else{
                               psi.new <- psi.b*(psi + mc.mle)
                               mc.var <- NA
                           }
                       })
        Psi[[i]] <- psi.new
        MC.datas[[i]] <- mc.data
        MC.MLEs[[i]] <- mc.mle
        MC.Hess[[i]] <- mc.Hess
        MC.Vars[[i]] <- mc.var

        psi <- psi.new
        flag.converge1 <- (mcl[1] < 1) & (mle.flag)
        flag.converge <- flag.converge1 & (n.s > n.s0)
        t.total <- t.total + tt[3]
        if(con$verbose){
            print(paste0(c("i = ", "n.s = ", "mc.mle = ", "converge = ", "time = "),
                         c(i, n.s, psi.new[1], flag.converge, t.total)))
        }
        n.s <- ifelse(flag.converge1 & (n.s == n.s0),  n.s*n.increase, n.s)
        mc.con$N.Zy <- n.s * mc.con$thin + mc.con$burns
        i <- i+1
    }
    if(con$trace.all){
        ans <- list(Psi=Psi, MC.datas=MC.datas, MC.MLEs=MC.MLEs, MC.Hess=MC.Hess,
                    MC.Vars=MC.Vars, data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "OptimMCL"
        return(ans)
    }else{
        ans <- list(Psi=Psi[[i-1]], MC.datas=MC.datas[[i-1]], MC.MLEs=MC.MLEs[[i-1]],
                    MC.Hess=MC.Hess[[i-1]], MC.Vars=MC.Vars[[i-1]],
                    data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "OptimMCL"
        return(ans)
    }
}

#### wrapper for direct and latent models
OptimMCL <- function(data, psi0, family, control = list(), mc.control = list()){
    if(family == "gauss"){
        OptimMCL.lm(data = data, psi0 = psi0, control = control)
    }else{
        OptimMCL.glm(data = data, psi0 = psi0, family = family, control = control, mc.control = mc.control)
    }
}
