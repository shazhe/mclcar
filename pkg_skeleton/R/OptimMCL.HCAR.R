#### iterative maximization of the Monte Carlo likelihood for the HCAR
OptimMCL.HCAR <- function(data, psi0,  control = list()){
    ## general algorithm default control
    con <- list(n.iter = 20, n.samples = 500, s.increase = 2, psi.ab = c(1,0.5),
                psi.range = c(-0.304, 0.0463, -0.3, 0.1737, 10, 10), mc.var = FALSE,
                trace.all = TRUE, verbose = TRUE, print = 1)
    con[names(control)] <- control

    ## Initializing
    n.s <- n.s0 <- con$n.samples
    n.increase <- con$s.increase
    n.iter = con$n.iter
    psi.a <- con$psi.ab[1]              # constant a in updating criterion 3.6
    psi.b0 <- con$psi.ab[2]              # constant b in updating criterion 3.6
    psi <- psi0


  if (is.environment(data)) {
    W <- get("W", envir=data)
    bb <- get("bb", envir=data)
  } else {
    W <- data$W
    bb <- data$bb
  }
    if(is.null(W)){
        rbb <- 1/range(bb)
        psi.range <- con$psi.range # range of rhos and range of scale change of sigmas
        ## Constraints on parameters
        A <- matrix(0, nrow = 6, ncol = length(psi0))
        A[,1] <- c(1, -1, 0, 0, 0, 0) # for lambda
        A[,2] <- c(0, 0, 1, -1, 0, 0) # for sigma_e
        A[,3] <- c(0, 0, 0, 0, 1, -1) # for sigma_u
        B <- c(-psi.range[3], psi.range[4],
               -psi[2]/psi.range[5], psi[2]*psi.range[5], -psi[3]/psi.range[6], psi[3]*psi.range[6])
    }else{
         psi.range <- con$psi.range # range of rhos and range of scale change of sigmas
        ## Constraints on parameters
        A <- matrix(0, nrow = 8, ncol = length(psi0))
        A[,1] <- c(1, -1, 0, 0, 0, 0, 0, 0) # for rho
        A[,2] <- c(0, 0, 1, -1, 0, 0, 0, 0) # for lambda
        A[,3] <- c(0, 0, 0, 0, 1, -1, 0, 0) # for sigma_e
        A[,4] <- c(0, 0, 0, 0, 0, 0, 1, -1) # for sigma_u
        B <- c(-psi.range[1], psi.range[2], -psi.range[3], psi.range[4],
               -psi[3]/psi.range[5], psi[3]*psi.range[5], -psi[4]/psi.range[6], psi[4]*psi.range[6])
    }

    Psi <- list()
    MC.datas <- list()
    MC.MLEs <- list()
    MC.Hess <- list()
    MC.Vars <- list()
    flag.converge <- FALSE
    i = 1
    t.total <- 0
    rtol = 1
    pscale <- rep(1, length(psi))
    pscale[2:3] <- 10
    pscale[5] <- 100

    while(i <= n.iter & !flag.converge){
        tt <- system.time({mc.data <- sim.HCAR(psi, data, n.samples = n.s)
              opt.res <- maxBFGS(fn = mcl.HCAR, grad = mcl.grad.HCAR,
                                 hess = mcl.Hessian.HCAR,
                               mcdata = mc.data, data = data,
                               start=as.numeric(psi), print.level=con$print,
                               parscale = pscale,
                               constraints=list(ineqA=A, ineqB=B), control = list(reltol = rtol))
            rtol <- max(opt.res$maximum/sqrt(n.s), sqrt(.Machine$double.eps))

            mc.mle <- opt.res$estimate
            mc.Hess <- opt.res$hessian
            mcl <- mcl.HCAR(mc.mle, mcdata = mc.data, data, Evar = TRUE)
            mcl.var <- mcl[2]
            mle.flag <- ifelse(is.na(mcl.var), FALSE, round(mcl.var*n.s) <= psi.a)
            if(mle.flag){
                psi.new <- mc.mle
                if(con$mc.var){
                    mc.var <- vmle.HCAR(opt.res, mcdata = mc.data, data = data)
                } else{
                    mc.var <-NA
                               }
            }else{
                #psi.b <- max(1 - 1/abs(mcl[1]), psi.b0)
                psi.new <- (1-psi.b0)*psi + psi.b0*mc.mle
                mc.var <- NA
            }
        })
        Psi[[i]] <- psi.new
        MC.datas[[i]] <- mc.data
        MC.MLEs[[i]] <- mc.mle
        MC.Hess[[i]] <- mc.Hess
        MC.Vars[[i]] <- mc.var

        psi <- psi.new
        flag.converge1 <- (mcl[1] < 3) & (mle.flag)
        flag.converge <- flag.converge1 & (n.s > n.s0)
        t.total <- t.total + tt[3]
        if(con$verbose){
            if(is.null(W)){
            message(paste0(c("i = ", "n.s = ", "lambda = ", "sigma.e = ",
                             "sigma.u = ", "beta0 = ", "converge = ", "time = "),
                           c(i, n.s, psi.new[1:4], flag.converge, t.total),
                           c(rep(c(",\t", ",\t", "\n"), 2), c(",\t", "\n"))))
        }else{
             message(paste0(c("i = ", "n.s = ", "rho = ", "lambda = ", "sigma.e = ",
                             "sigma.u = ", "beta0 = ", "converge = ", "time = "),
                           c(i, n.s, psi.new[1:5], flag.converge, t.total),
                           c(rep(c(",\t", "\n"), 3), "\n", c(",\t", "\n"))))
        }
        }
        n.s <- ifelse(flag.converge1 & (n.s == n.s0),  n.s*n.increase, n.s)
        i <- i+1
    }
    if(con$trace.all){
        ans <- list(Psi=Psi, MC.datas=MC.datas, MC.MLEs=MC.MLEs, MC.Hess=MC.Hess,
                    MC.Vars=MC.Vars, data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "OptimMCL.HCAR"
        return(ans)
    }else{
        ans <- list(Psi=Psi[[i-1]], MC.datas=MC.datas[[i-1]], MC.MLEs=MC.MLEs[[i-1]],
                    MC.Hess=MC.Hess[[i-1]], MC.Vars=MC.Vars[[i-1]],
                    data=data, N.iter = i-1, total.time = t.total,
                    convergence = flag.converge, mcsamples = c(n.s0, n.s))
        attr(ans, "class") <- "OptimMCL.HCAR"
        return(ans)
    }
}
